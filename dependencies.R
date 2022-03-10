# No Remotes ----
# Attachments ----
to_install <- c("config", "golem", "processx", "shiny","shinydashboard","shinyWidgets","plotly","janitor","tidyverse","DT","openxlsx","WriteXLS","waiter","ggpubr","reshape2","fastmatch","ggthemes","xml2")
  for (i in to_install) {
    message(paste("looking for ", i))
    if (!requireNamespace(i)) {
      message(paste("     installing", i))
      install.packages(i)
    }
  }
