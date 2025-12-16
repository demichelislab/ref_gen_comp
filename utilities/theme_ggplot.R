theme_publication <- function(){ 
  font <- "Helvetica"   #assign font family up front
  
  theme_bw(base_size = 15) %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      
      #text = element_text(size = 15),
      
      #text elements
      plot.title = element_text(             #title
        family = font,            #set font family
        size = 27,                #set font size
        face = 'bold',            #bold typeface
        hjust = 0.5,                #left align
        vjust = 1),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
        family = font,            #font family
        size = 25,
        hjust = 0.5,
        vjust = 1),               #font size
      
      plot.caption = element_text(           #caption
        family = font,            #font family
        size = 20,                 #font size
        hjust = 1),               #right align
      
      # panel.background = element_rect(fill = "white"), # necessary to avoid drawing panel outline
      # plot.background = element_rect(fill = "white",
      #                                colour = NA), # necessary to avoid drawing plot outline
      # legend.background = element_rect(fill = "white", colour = NA),
      # legend.box.background = element_rect(fill = "white", colour = NA),
      # legend.key = element_rect(fill = "white", colour = NA)
      # # axis.title = element_text(             #axis titles
      #   family = font,            #font family
      #   size = 24),               #font size
      # 
      # axis.text = element_text(              #axis text
      #   family = font,            #axis famuly
      #   size = 18),                #font size
      # 
      # strip.text = element_text(
      #   family = font, 
      #   size = 18,
      #   margin(l=1, b =1, r=1, t=1)),
      # 
      # legend.title = element_text(              #axis text
      #   family = font,            #axis famuly
      #   size = 22),
      # 
      # legend.text = element_text(              #axis text
      #   family = font,            #axis famuly
      #   size = 20)
      #,
      
     # legend.position = "bottom"
    )
}
