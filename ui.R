#load libraries
library(DT)
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinycssloaders)

#Hack to get rid of  dashboard button
mydashboardHeader <- function(..., title = NULL, disable = FALSE,title.navbar=NULL, .list = NULL) {
  items <- c(list(...), .list)
  tags$header(class = "main-header",
  style = if (disable) "display: none;",
  span(class = "logo", title),
  tags$nav(class = "navbar navbar-static-top", role = "navigation",
  span(shiny::icon("bars"), style = "display:none;"), 
  title.navbar,
  div(class = "navbar-custom-menu",
  tags$ul(class = "nav navbar-nav",items)
  )))
}
#===================================================================================================================================================================================== 
#Ok Defining UI page 
ui <- dashboardPage(
#define header title           
  mydashboardHeader( title=tags$p("ID-TaxER",style="font-size: 40px; font-family:Calibri; background-color: #337ab7;")),
#width of dashboard sidebar              
  dashboardSidebar(width=230, 
#logo image matching width to sidebar width
  img(height=231.141,width=230,src="logos.png"),
#more info button, again keeping size in mind
  actionButton("more_info_button", "More Information",   style="color: #fff; background-color: #337ab7;
                  border-color: #2e6da4",width=200),
#github button
  actionButton("github_button", "  GitHub Page  ",onclick ="window.open('https://github.com/brijon/ID-TaxER-flat-files', '_blank')",
                 style="color:#fff; background-color: #337ab7; border-color: #2e6da4",width=200)
  ),               
  dashboardBody(
    tags$head(
    tags$style(HTML(
        '.myClass { 
        font-size: 25px;
        line-height: 50px;
        text-align: center;
        font-family: "Calibri";
        padding: 0 130px;
        overflow: hidden;
        color: white;
        }'))
    ),
    tags$script(HTML('
#adds title to top of page
    $(document).ready(function() {
    $("header").find("nav").append(\'<span class="myClass"> Identification of Taxa & Environment Responses </span>\');
    })
    ')),
#Hidden section about app
    useShinyjs(),
#tags$h is heading, tags$p is paragraph etc , tags$b is bold etc                
    shinyjs::hidden(div(id="more_info",tags$h4(tags$b("Summary"),style="color:#000080;font-family:Calibri;"),
    tags$p( "ID-TaxER provides an interface to explore potential soil habitat preferences of bacterial taxa derived from 16S rRNA gene sequencing. Query sequences are blasted against a database of representative sequences of 97% OTUs obtained from a large soil survey conducted across Britain (the Countryside Survey). Each sequence in the database is linked to an additional trait matrix containing taxonomic assignments as well as environmentally derived information about that OTU (e.g pH or habitat preference). Results are displayed as an interactive table of hits with percentage match to  a CS sequence, and associated taxonomy (greengenes). Upon selecting a hit, a plot of model fit to soil parameters is displayed indicating for example the pH optima of that taxon, as well as habitat preferences and spatial distribution (currently Britain only).  
		",style="color:#000080;font-family:Calibri;"),
    tags$h4(tags$b("Limitations"),style="color:#000080;font-family:Calibri;"),
		tags$p("The database encompasses the V3-V4 region of the 16SrRNA, amplified with 341f/806r primers. Queries which do not cover this region will obviously give incorrect results, and additionally taxa poorly amplified with these primers will be under represented. Importantly this tool is based on homology mapping to a short portion of the conserved 16S rRNA gene, and so all the usual limitations apply regarding accuracy of taxonomic (and habitat preference) assignment.  It is therefore for research purposes only.",style="color:#000080;font-family:Calibri;"),
		tags$h4(tags$b("Ongoing work"),style="color:#000080;font-family:Calibri;"),
		tags$p("A",tags$a(href="https://github.com/brijon/ID_TaxER-Custom-Database-for-DADA2","github",target="_blank")," page has been set up for batch sequence querying using Dada2. We are exploring options to allow users to upload their own ecological trait information to the trait matrix (e.g if a sequence with high homology to a CS sequence comes up in user experiments as an indicator of warming, drought, plant species X etc then it would be useful to capture this information in the trait matrix).
		We also have similar ITS and 18S datasets which could be developed in a similar portal if enough interest.",style="color:#000080;font-family:Calibri;"),
		tags$h4(tags$b("Contact:"),style="color:#000080;font-family:Calibri;"),
		tags$p("Briony Jones", tags$b("(brijon@ceh.ac.uk)"),style="color:#000080;font-family:Calibri;"),
    tags$p("Rob Griffiths",tags$b(" (rig@ceh.ac.uk)
    "),style="color:#000080;font-family:Calibri;"))),
    br(),
#enter sequence box  place holder initially empty    
    textInput(inputId = "mysequence",label="Please enter a sequence",
    value="",width = 10000, placeholder = ''),
    br(),
#various buttons               
    actionButton("blast", "Blast",style="color: #fff; background-color:#B20000"),
    actionButton("resetSequence", "Clear Input",style="color: #fff; background-color:#228B22"),
    actionButton("exampleSequence", "Example Sequence",style="color: #fff; background-color:#98AFC7"),
    br(),
#area where any warning messages appear                
    span(textOutput("Warning"),style="color:red;font-size:17px  "),
    br(),
#hidden results section            
    shinyjs::hidden(div(id="Results", hr(),
#===================================================================================================================================================================================== 
#hits table
    fluidRow(column(width=8,hr(),HTML('<center><h4>Top Hits</h4></center>'),box(DT::dataTableOutput("blastout"),width=500,height=578)),
#===================================================================================================================================================================================== 
#define collumn ie portion of interface to the right for plots etc                
    column(width=4,hr(),
    tabsetPanel(id="plotTabset",
#===================================================================================================================================================================================== 
#first plot HOF                
    tabPanel(title="pH Model",box(id="plotbox",plotOutput('modelplot'), width=350,height=425)
#end of plot tab panel
    ),
#===================================================================================================================================================================================== 
#LOESS plot              
    tabPanel(title="pH LOESS",box(id="plotbox2",plotOutput('loessplot'), width=350,height=425)
#end of plot tab panel
    ),
#===================================================================================================================================================================================== 
#Map                
    tabPanel(title="GB Map",box(id="Map_box",withSpinner(plotOutput('map'),type=7),width=350,height=425)
#	end of plot tab panel
	  ),                  
#===================================================================================================================================================================================== 
#AVC plot
    tabPanel(title="Habitats",box(id="AVC_plotbox",plotOutput('AVC_box_plot'),width=350,height=425)
    # end of plot tab panel
    )
#===================================================================================================================================================================================== 
#indicators text box
    ,fluidRow( HTML('<center><h4>Additional Information (User Submitted)</h4></center>'),box(
 #want to make scrollable box               
    tags$style(HTML("
      #Indicators {
        height:64px;
         overflow-y:scroll
        }")
    ),
    width=50,
    height=80,  
    htmlOutput('Indicators')
    ))
#===================================================================================================================================================================================== 
#end of tabset(plots and indicators)
    )
#end of collumn
    )  
#===================================================================================================================================================================================== 
#blast output and taxonomy
  ,hr()
  ,fluidRow(column(width=12,HTML('<h4><center>Blast Output and Taxonomy</center></h4>'),box(DT::dataTableOutput("BlastResults"),DT::dataTableOutput("OTU.Taxon")
  ,width=1000)
  )
#end of fluid row
  )
#===================================================================================================================================================================================== 
#end of hidden results divider
  )
#end of dashboard body
  )
#end of dashboard page
)))
