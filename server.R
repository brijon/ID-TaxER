#Load libraries
library(shiny)
library(shinyjs)
library(DT)
library(RPostgreSQL)
library(pool)
library(eHOF)
library(maptools)
library(vegan)
library(ggplot2)
library(gstat)
library(RColorBrewer)
library(reshape)
#===================================================================================================================================================================================== 
#define server
server <- function(input, output,session) {
	
#connect to postgres database
  con <- dbPool(
 	  drv = RPostgreSQL::PostgreSQL(max.con=140),
 	  dbname = "molecular",
 	  host = "shiny-prod.nerc-lancaster.ac.uk",
 	  port='5432',
 	  user =getOption("trait_matrix_userid"),
	  password =getOption("trait_matrix_password")

  )
#===================================================================================================================================================================================== 
#shiny JS commands
#get info about app , shiny js allows "more information" to appear and dissapear when clicking more information         	
  onclick("more_info_button",toggle(id="more_info",anim=TRUE))
#multiple show functions so specify shinyjs version 
  onclick("blast",shinyjs::show(id="Results",anim=TRUE))
#===================================================================================================================================================================================== 
#get all sample environment information going to use this for habitat box plots
  SQL_command=paste("select * from env_attributes.env_attributes_ur;")
#env consists of sample avc_code (habitat code), avc (habitat description) and pH
  env <- dbGetQuery(con, SQL_command)
#make sample names row names
  rownames(env)<-env[,1]
#correct values with "NA" instead of NA  
  env[env=="NA"]<-NA
#for purpose of loading serialised objects 
  dbGetQuery(con, "set standard_conforming_strings to 'on'")
#get uk map outline to plot uk mapping objects  
  SQL_command=paste("select plot_object from plotting_tools.map_tools where description= 'map_outline';")
  uk.line <-unserialize(postgresqlUnescapeBytea( dbGetQuery(con, SQL_command))) 
#===================================================================================================================================================================================== 
#define function that runs blast (sequence allignment) ,arguments are blast 16S database and query sequence
  make.comparison <- function(db, query ){
#check query isnt empty
    if (query!=""){    
#blast_command for alligning sequences returmns top 20 hits
      cmd <- paste('echo -e ">Name\n',query,'"', '|/home/brijon/ncbi-blast-2.7.1+/bin/blastn' ,'-db',
'/home/brijon/repseqs_ID_TaxER_db/repseqs_ID_TaxER','-num_alignments',20,'-evalue',0.001, '-outfmt', 7)
   
#run system command and capture output
      blast_capture<- system(paste("/bin/bash -c", shQuote(cmd)),intern=TRUE)
#check there are hits
      if(blast_capture[4]!="# 0 hits found"){
#this variable will be used later to identify that output should be displayed (as hits have been returned from blast command)      
        output_switch <-"on"
#remove first and last lines of output (descriptive not hits)
        blast_capture<-blast_capture[-length(blast_capture)]
        blast_capture<-blast_capture[6:25]
#some dont have 20 hits! so nas appear get rid of them 
        blast_capture<-blast_capture[!is.na(blast_capture)]  
#===================================================================================================================================================================================== 
#make empty arrays/dataframes
#this will list all OTU's the sequence has hit to 
        OTUlist<-c()
#this table will become main results table on ID_TaxER with subject_ID, kingdom, phylum, class, order, family, genus and species     
#empty matrix
        m <- matrix(nrow=0, ncol = 9)
#convert to dataframe
        ph.model<-data.frame(m)  

#tax dataframe this will display underneath main table 
        m<-matrix(nrow=0, ncol=8)
        tax<-data.frame(m)
        
#this variable will be used for full blast output that will also display underneath main table         
        all_blst_output<-c()

#abund dataframe will be used to store microbial relative abundance for plotting habitat boxplots      
#get col names for abundance data (sample names)
        myQuery <- "SELECT * FROM information_schema.columns WHERE table_schema = 'otu_abund' AND table_name = 'otu_abund_ur'"
        samp_colnames <- dbGetQuery(con, myQuery)[,4]
 
#make empty matrix with sample names as colnames     
        m<-matrix(nrow=0,ncol=1007)
        abund<-data.frame(m)
        colnames(abund)<-samp_colnames
#===================================================================================================================================================================================== 
#fill arrays and df
        for (blast_hit in blast_capture){
#split blast fields by tab 
          blast_hit_fields=strsplit(blast_hit,"\t")[[1]]
          OTU=blast_hit_fields[2]
# look up OTU in postgres ph_model_features table this gives us OTU hit, hof model,optimum1,optimum2,PpH reclass1, pH reclass2, abundance_rank and occupancy_proportion
          SQL_command=paste("select * from model_attributes.ph_model_features_ur WHERE hit='",toString(OTU),"';",sep="")
          ph.model.row <- dbGetQuery(con, SQL_command)
#ok lets add blast percentage Identity to this         
          ph.model.row$Blast_Percentage_Identity<-blast_hit_fields[3]
#attach to empty ph.model df and put blast identity first 
          ph.model<-rbind(ph.model,ph.model.row[,c(9,1:8)])
#remove duplicates if different areas of overlap of same sequence     
          ph.model<- subset(ph.model, !duplicated(ph.model$hit))
#make OTU hit rownames         
          row.names(ph.model)<-ph.model$hit
#Now look up OTU in postgres to get tax info	     
          SQL_command=paste("SELECT * FROM otu_attributes.taxonomy_ur WHERE hit='",toString(OTU),"';",sep="")
#should now have row with hit, kingdom, phylum, class, order, family, genus and species
          tax.row<-dbGetQuery(con,SQL_command)
#attach to tax df
          tax<-rbind(tax,tax.row)
          tax<- subset(tax, !duplicated(tax$hit))
#look up otu in postgres to get abundance info    
          SQL_command=paste("SELECT * FROM otu_abund.otu_abund_ur WHERE hit='",toString(OTU),"';",sep="")
          abund.row<-dbGetQuery(con,SQL_command)
#attach to abund df      
          abund<-rbind(abund,abund.row)
          abund<- subset(abund, !duplicated(abund$hit))
#attach to all blast_output     
          all_blst_output<-c(all_blst_output,strsplit(blast_hit,"\t"))
          OTUlist=c(OTUlist,OTU)
    
        }
#===================================================================================================================================================================================== 
#remove any replicates
        OTUlist<-unique(OTUlist)
#make blast output into nice df for the GUI 
        all_blst_output_df<-t(as.data.frame(all_blst_output))[,2:12]
#make OTU's rownames
        row.names(all_blst_output_df)<-all_blst_output_df[,1]
#make colnames
        colnames(all_blst_output_df)<-c('subject id',' % identity', 'alignment length', 'mismatches', 'gap opens', 'q. start', 'q. end', 's. start', 's. end', 'evalue', 'bit score')
        all_blst_output_df<-as.data.frame(all_blst_output_df[OTUlist,])
#remove duplicates from all blast output... could happen if different areas of overlap of same sequence
#assign rownames to tax df	
        row.names(tax)<-tax$hit
#assign colnames to taxonomy df
        colnames(tax)=c("subject id","kingdom","phylum","class","order","family","genus","species")
        tax<-as.data.frame(tax)
#add colnames to ph.model dataframe      
        colnames(ph.model)<-c("Blast Percentage Identity","CS OTU hit","pH HOF model","pH Optimum 1","pH Optimum 2","Model description","pH Class","Abundance rank","Occupancy")
        ph.model<-as.data.frame(ph.model)
        abund<-t(abund)
#make first row of abund colnames      
        colnames(abund)<-abund[1,]
        abund<-abund[-1,]
#return all tables      
        return(list(ph.model,all_blst_output_df,tax,abund,output_switch))
      }
      else{
#if query has no hits, switch variable is given value off to identify that output should not be displayed         
        output_switch <-"off"
        return(output_switch)       
      }
    }
    else{
#if query empty,switch variable is given value off to identify that output should not be displayed     
      output_switch <-"off"
      return(output_switch)

    }
  }
#===================================================================================================================================================================================== 
 #make blast run when you click blast button using eventReactive Function
 run_sequence<-eventReactive(input$blast,{make.comparison("/home/brijon/repseqs_db", input$mysequence) })
#===================================================================================================================================================================================== 
#define main table 
#selection mode is one at a time (dont want to show several plots at once etc), automatically select top row otherwise other outputs will be empty when page loads which is just ugly :)
  output$blastout<-DT::renderDataTable({
    all_mod_output<-run_sequence()
#if switch veriable is on and has hits
   if(all_mod_output[[length(all_mod_output)]]=="on"){
     ph_mod_output<-all_mod_output[[1]]
   }
   }, selection = list(mode='single',selected=1),options=list(scrollX=TRUE,pageLength=7,dom='tp'),rownames=FALSE, colnames = c("Blast Percentage Identity","CS OTU hit","pH HOF model","pH Optimum 1","pH Optimum 2","Model description","pH Class","Abundance rank","Occupancy")
  ) 
#===================================================================================================================================================================================== 
#HOF plot code
  output$modelplot = renderPlot({
    all_mod_output<-run_sequence()
#if switch variable on and ther are hists
    if(all_mod_output[[length(all_mod_output)]]=="on"){
      s=input$blastout_rows_selected
#get selected rows in df
      if (length(s)){
        OTU=all_mod_output[[2]][s,1]
        par(mar = c(4, 4, 1, 4))
        dbGetQuery(con, "set standard_conforming_strings to 'on'")
#get HOF model object        
        ph.mod.obj=dbGetQuery(con, paste("SELECT mod_object FROM models.ph_ur WHERE hit='",toString(OTU),"';",sep="")) 
#unserialize 
        ph.mod.obj<-postgresqlUnescapeBytea(ph.mod.obj)
        ph.mod.obj<- unserialize(ph.mod.obj)
      
#own custom plot with relative abundance as modifying plot hof object proved difficult 
        mod_choice<-Para(ph.mod.obj)$model
#get model fit         
        fitted=ph.mod.obj$models[mod_choice][[1]][9] 
        mod_stats=as.data.frame(cbind(ph.mod.obj$y,ph.mod.obj$x,fitted$fitted))
        mod_stats=mod_stats[order(mod_stats$V2),]
        mod_col=c("black","red","#2BF33F","#3554EE","#895A3B")
        names(mod_col)=c("I","II","III","IV","V")
#plot          
        plot(mod_stats$V2,mod_stats$V1,ylim=c(0,quantile(mod_stats$V1,0.999)),cex=0.7,pch=19,xlab="pH",ylab=paste("Number of reads ","(",toString(OTU),")",sep=""))
#add model fit         
        lines(mod_stats$V2,mod_stats$V3,ylim=c(0,quantile(mod_stats$V1,0.999)),lwd=2,col=mod_col[mod_choice])
      }
    }
    else{
#if output switch off     
      plot.new()
    }
#end of render plot   
  })
#===================================================================================================================================================================================== 
#LOESS plot code 
  output$loessplot = renderPlot({
    all_mod_output<-run_sequence()
#if output switch on 
    if(all_mod_output[[length(all_mod_output)]]=="on"){
#selected rows 
      s=input$blastout_rows_selected
      if(length(s)){
        par(mar = c(4, 4, 1, 4)) 
        abund <-all_mod_output[[4]]
#get relevant rows from abundance table          
        OTU.spc<-as.data.frame(as.numeric(abund[,s]),row.names=row.names(abund))
        OTU.spc_pH<-merge(OTU.spc,env,by=0)
#drop unecessary collumns
        OTU.spc_pH<-OTU.spc_pH[,-c(3:5)]
#order df by ph
        OTU.spc_pH<-OTU.spc_pH[order(OTU.spc_pH[,3]),]
#calculate loess         
        loessMod50 <- loess(OTU.spc_pH[,2]~OTU.spc_pH[,3], span=0.5)
        smoothed50 <- predict(loessMod50,se=TRUE) 
#plot points         
        plot(OTU.spc_pH[,3],OTU.spc_pH[,2],ylim=c(min(OTU.spc_pH[,2]),quantile(OTU.spc_pH[,2],0.999)),cex=0.7,pch=19,xlab="pH",ylab=paste("Relative Abundance (",colnames(abund)[s],")",sep=""),col="#878273")
#plot loess         
        lines(smoothed50$fit, x=OTU.spc_pH[,3],ylim=c(min(OTU.spc_pH[,2]),quantile(OTU.spc_pH[,2],0.999)), col="black",lwd=2)
      }
    }     
  })     
#===================================================================================================================================================================================== 
#user submitted section (indicators, e.g additional information about sequence)
  output$Indicators=renderText({
    all_mod_output<-run_sequence()
#if switch variable is set to on          
    if(all_mod_output[[length(all_mod_output)]]=="on"){
      s=input$blastout_rows_selected
#get selected rows in df
      if (length(s)){
        OTU=all_mod_output[[2]][s,1]
        SQL_command=paste("select * from otu_attributes.user_submitted_ur WHERE hit='",toString(OTU),"';",sep="")
        indicator.tab<- dbGetQuery(con, SQL_command)
#if there are indicators for this otu                
          if (nrow(indicator.tab)>=1){
#empty variable to add indicators to              
            indicators=""
#loop over difference datasets in reference collumn
            for (study in unique(indicator.tab$reference)){
#get indicator subsets  for each dataset                  
              study_subset=indicator.tab[which(indicator.tab$reference==study),]
#sort by identity and pvalue
              study_subset=study_subset[order(-study_subset$identity,study_subset$pvalue),]
              indicator=paste(study_subset[1,4],"<b>  p value: </b>",study_subset[1,5],"\n <b> Source: </b>",study_subset[1,6],sep="") 
#paste with line breaks
              indicators=paste(indicator,indicators,sep='<br/>')
            }
            indicators=paste("<center> <b>",OTU,"</center> </b> <center>",indicators,"</center>",sep="")
   		      HTML(indicators)
          }else{
            noindicators= paste("<center> <b>",OTU," </b> </center> \n <center> No user submitted information.</center>",sep="")
            HTML(noindicators)}
      }
    }
  })

#===================================================================================================================================================================================== 
#map
  output$map=renderPlot({
#if switch on     
    all_mod_output<-run_sequence()
    if (all_mod_output[[length(all_mod_output)]]=="on"){
      s=input$blastout_rows_selected
      if(length(s)){
	      	OTU=all_mod_output[[2]][s,1]
          par(mar = c(4, 4, 1, 4))
          dbGetQuery(con, "set standard_conforming_strings to 'on'")
#unserialise          
          map.obj=dbGetQuery(con, paste("SELECT map_object FROM otu_attributes.maps_ur WHERE hit='",toString(OTU),"';",sep=""))
          raw_retrieved_map<-postgresqlUnescapeBytea(map.obj) 
       		object_retrieved_map<-unserialize(raw_retrieved_map)
#plot       		
          spplot(object_retrieved_map[[1]]["var1.pred"],at=unlist(object_retrieved_map[2:11]),xlab=toString(OTU),ylab.pos=c(5,10,100),sp.layout=list("sp.lines",uk.line,lwd=2,col="black"))
      }
    }
  })
#===================================================================================================================================================================================== 
  output$AVC_box_plot = renderPlot({    
    all_mod_output<-run_sequence()
#if swithch var on 
    if(all_mod_output[[length(all_mod_output)]]=="on"){
#selected rows 
      s=input$blastout_rows_selected
      if(length(s)){
        abund <-all_mod_output[[4]]
#get relevant abund info         
         OTU.spc<-as.data.frame(as.numeric(abund[,s]),row.names=row.names(abund))
         numb_otus<-specnumber(OTU.spc)
         if (sum(numb_otus)>150){
          ord<-reorder(env$avc,X=as.numeric(env$avc_code),FUN=mean)
          par(mar=c(11,5.8,2,0.4)+0.1)
          boxplot(OTU.spc[,1]~ord,las=2,ylab="",cex=0.5,outline=FALSE)
          title(ylab=paste("Relative Abundance (",colnames(abund)[s],")",sep=""), line=4.5, cex.lab=0.9)
         }else{
          par(mar = c(0,0,0,0))
          empty_plot=plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
          empty_plot+ text(x = 0.8, y = 0.8, paste("Not enough data to generate box plot"),
          cex = 1, col = "red", family="sans", font=1, adj=1)
        }
      }
    }
    else{
      plot.new()
    }
  })
#===================================================================================================================================================================================== 
#taxonomy
  output$OTU.Taxon<-DT::renderDataTable({
    all_mod_output<-run_sequence()
    if(all_mod_output[[length(all_mod_output)]]=="on"){
      s=input$blastout_rows_selected
      if (length(s)) {
        as.data.frame(all_mod_output[[3]])[s,]
      }
    }
#no selection and just tables no extras like search(dom=t)
  },selection='none',options=list(dom='t'),rownames=FALSE)
#===================================================================================================================================================================================== 
#blast output
  output$BlastResults<-DT::renderDataTable({
#dat_raw
    all_mod_output<-run_sequence()
    if(all_mod_output[[length(all_mod_output)]]=="on"){
      s=input$blastout_rows_selected
      if (length(s)){
        all_mod_output[[2]][s,]
#end of second if statement    
      }
   #end of first if statement
    }
     #end of render datatable 
  },selection='none',options=list(dom='t'),rownames=FALSE)
#displays if switch set to off  
  output$Warning <- renderText({ 
    all_mod_output<-run_sequence()
    if(all_mod_output[[length(all_mod_output)]]!="on"){
      paste("No Hits Found!")
    }
     
  })

#if reset button is pressed, query is reset   
  observeEvent(input$resetSequence, {
    reset("mysequence")
  })

#if example sequence is selected, example sequence entered into query box
  observeEvent(input$exampleSequence, {
    updateTextInput(session,"mysequence",value="ACAGAGGTCTCAAGCGTTGTTCGGATTCATTGGGCGTAAAGGGTGCGTAGGTGGTGATGCAAGTCTGGTGTGAAATCTCGGGGCTCAACTCCGAAATTGCACCGGATACTGCGTGACTCGAGGACTGTAGAGGAGATCGGAATTCACGGTGTAGCAGTGAAATGCGTAGATATCGTGAGGAAGACCAGTTGCGAAGGCGGATCTCTGGGCAGTTCCTGACACTGAGGCACGAAGGCCAGGGGAGCAAACGGG")
  })  
}



     
