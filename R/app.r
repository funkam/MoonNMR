#' The application User-Interface and server
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny shinydashboard shinyWidgets plotly ggplot2 waiter ggpubr reshape2 fastmatch xml2 dplyr tidyr ggfortify
#' @importFrom WriteXLS WriteXLS
#' @importFrom openxlsx read.xlsx
#' @importFrom stringr str_c str_detect
#' @importFrom janitor clean_names remove_constant
#' @importFrom stats prcomp
#' @importFrom utils read.csv read.csv2 unzip write.csv
#' @importFrom grDevices dev.off pdf
#' 
#' @noRd


startProjects<-c()
#Script group_plotter----
plotter_grouped<-function(data,group,columns,stats) {
  data[data==0] <- NA
  columns<-c(group,columns)
  data<-select(data,all_of(columns))
  metabolites<-names(data)
  varx<-metabolites[1]
  
  
  if(stats=="none"){
    pltg<- lapply(seq_len(ncol(data)), FUN = function(x) {
      ggplot(data,aes(x = .data[[varx]], y = data[ , x],inheriet.aes=FALSE,fill=.data[[varx]],group=.data[[varx]])) +
        geom_boxplot()+
        theme_classic()+
        theme(legend.position="none",axis.text.x=element_text(angle=45,hjust=1))+
        labs(x="", y="", title = colnames(data)[x])
    })
  }
  if(stats=="ttest"){
    pltg<- lapply(seq_len(ncol(data)), FUN = function(x) {
      ggplot(data,aes(x = .data[[varx]], y = data[ , x],inheriet.aes=FALSE,fill=.data[[varx]],group=.data[[varx]])) +
        geom_boxplot()+
        theme_classic()+
        theme(legend.position="none",axis.text.x=element_text(angle=45,hjust=1))+
        labs(x="", y="", title = colnames(data)[x])+
        stat_compare_means(method="t.test",aes(label=..p.signif..),hide.ns = TRUE)
    })
  }
  if(stats=="wilcox"){
    pltg<- lapply(seq_len(ncol(data)), FUN = function(x) {
      ggplot(data,aes(x = .data[[varx]], y = data[ , x],inheriet.aes=FALSE,fill=.data[[varx]],group=.data[[varx]])) +
        geom_boxplot()+
        theme_classic()+
        theme(legend.position="none",axis.text.x=element_text(angle=45,hjust=1))+
        labs(x="", y="", title = colnames(data)[x])+
        stat_compare_means(method="wilcox.text",aes(label=..p.signif..),hide.ns = TRUE)
      
    })
  }
  if(stats=="anova"){
    pltg<- lapply(seq_len(ncol(data)), FUN = function(x) {
      ggplot(data,aes(x = .data[[varx]], y = data[ , x],inheriet.aes=FALSE,fill=.data[[varx]],group=.data[[varx]])) +
        geom_boxplot()+
        theme_classic()+
        theme(legend.position="none",axis.text.x=element_text(angle=45,hjust=1))+
        labs(x="", y="", title = colnames(data)[x])+
        stat_compare_means(method="anova",aes(label=..p.signif..),hide.ns = TRUE)
    })
  }
  if(stats=="kruskal"){
    pltg<- lapply(seq_len(ncol(data)), FUN = function(x) {
      ggplot(data,aes(x = .data[[varx]], y = data[ , x],inheriet.aes=FALSE,fill=.data[[varx]],group=.data[[varx]])) +
        geom_boxplot()+
        theme_classic()+
        theme(legend.position="none",axis.text.x=element_text(angle=45,hjust=1))+
        labs(x="", y="", title = colnames(data)[x])+
        stat_compare_means(method="kruskal.text",aes(label=..p.signif..),hide.ns = TRUE)
    })
  }
  pltg<-pltg[-1]
  plot_combined<-ggarrange(plotlist=pltg)
  plot<-append(plot,list(plot_combined))
}






#Script xmler----
xmler<-function(data2,type){

  l<-data2

  if(type=="Metabolites"){

    #reading in the xml file
    all_measurements<-lapply(seq_len(length(data2)), FUN=function(x) {
      l<-data2[[x]]


      #numbner of parameters
      parameters_xml<-xml_find_all(l,".//VALUE")

      #Find the Sample name and clean it up
      Sample<-xml_find_all(l,".//SAMPLE")
      Sample<-xml_attr(Sample,"name")
      Sample<-sub("_e.*","",Sample)

      #create data frame
      xmlasdataframe<-as.data.frame(lapply(seq_len(length(parameters_xml)), FUN=function(x) {
        names<-xml_attrs(xml_child(xml_child(l, 4), x))[["name"]]
        conc<-xml_attrs(xml_child(xml_child(xml_child(l, 4), x), 1))[["conc"]]
        errConc<-xml_attrs(xml_child(xml_child(xml_child(l, 4), x), 2))[["errConc"]]
        rawConc<-xml_attrs(xml_child(xml_child(xml_child(l, 4), x), 2))[["rawConc"]]
        sigCorr<-xml_attrs(xml_child(xml_child(xml_child(l, 4), x), 2))[["sigCorr"]]
        lod<-xml_attrs(xml_child(xml_child(xml_child(l, 4), x), 1))[["lod"]]
        #loq<-xml_attrs(xml_child(xml_child(xml_child(l, 4), x), 1))[["loq"]]
        rbind(names,conc,errConc,rawConc,sigCorr,lod)
      }
      ))
      df1<-data.frame(t(xmlasdataframe[,-1]))
      rownames(df1)<-c()
      melty<-melt(df1,id.vars = "names")
      melty$names<-str_c(melty$names,"_",melty$variable)
      melty<-subset(melty,select=-c(variable))
      final_table<-melty[order(melty$name),]
      final_table<-pivot_wider(final_table,names_from="names",values_from="value")
      final_table<-cbind(Sample,final_table)
      if(x==1){
        thisisit<<-final_table[FALSE,]
      }
      thisisit<-rbind(thisisit,final_table)
      thisisactuallyit<-thisisit

    })
    test<-do.call(rbind.data.frame,all_measurements)
  }

  if(type=="Lipids"){

    all_measurements<-lapply(seq_len(length(data2)), FUN=function(x) {

      l<-data2[[x]]


      #numbner of parameters
      parameters_xml<-xml_find_all(l,".//VALUE")

      #Find the Sample name and clean it up
      Sample<-xml_find_all(l,".//SAMPLE")
      Sample<-xml_attr(Sample,"name")
      Sample<-sub("_e.*","",Sample)


      xmlasdataframe<-as.data.frame(lapply(seq_len(length(parameters_xml)), FUN=function(x) {
        names<-xml_attrs(xml_child(xml_child(l, 4), x))[["name"]]
        values<-xml_attrs(xml_child(xml_child(xml_child(l, 4), x), 1))[["value"]]
        unit<-xml_attrs(xml_child(xml_child(xml_child(l, 4), x), 1))[["unit"]]
        rbind(names,values,unit)
      }))
      df1<-data.frame(t(xmlasdataframe[,-1]))
      rownames(df1)<-c()
      df1$names<-str_c(df1$names,"_",df1$unit)
      df1<-subset(df1,select=-c(unit))
      #remove duplicated value inoriginal xml file
      df1<-df1[!duplicated(df1),]

      final_table<-df1[order(df1$name),]
      final_table<-distinct(pivot_wider(final_table,names_from="names",values_from="values"))
      final_table<-cbind(Sample,final_table)
      if(x==1){
        thisisit<<-final_table[FALSE,]
      }
      thisisit<-rbind(thisisit,final_table)
      thisisactuallyit<-thisisit

    })
    test<-do.call(rbind.data.frame,all_measurements)
  }
  test<-test


}


#Script Plotter_single-----
plotter_single<-function(data,group,column,grouped,stat){

  ata_original<-data
  if(length(column)>1){
    column<-column[1]
  }
  columns<-c(group,column)
  data<-select(data,all_of(columns))
  data[data==0] <- NA
  #if length ofcolumn vector >1 then show message....
    plot<-ggplot(data,aes(x = data[ ,group], y = data[ ,column],inherit.aes=FALSE,group= data[ ,group],fill= data[ ,group])) +
    geom_boxplot()+
    theme_classic()+
    theme(legend.position="none")+
    labs(x="", y="", title = column)
  plot<-plot
}
#SCRIPT ExperimentPicker----
experimentpicker<-function(solvent,size){

  if(solvent=="Urine"){
    if(size=="5mm"){
      EXPERIMENTE<-c("N PROF_URINE_NOESY","N PROF_URINE_JRES","N PROF_URINE_DAS_A","N PROF_URINE_DAS_E")
    }
    else if(size=="3mm"){
      EXPERIMENTE<-c("N PROF_URINE_NOESY_3mm","N PROF_URINE_JRES_3mm","N PROF_URINE_DAS_A","N PROF_URINE_DAS_E")
    }
  }
  if(solvent=="Urine_Neo"){
    if(size=="5mm"){
      EXPERIMENTE<-c("N PROF_URINE_NOESY","N PROF_URINE_JRES","N PROF_URINE_DAS_N")
    }
    else if(size=="3mm"){
      EXPERIMENTE<- c("N PROF_URINE_NOESY_3mm","N PROF_URINE_JRES_3mm","N PROF_URINE_DAS_N")
    }
  }
  if(solvent=="Plasma"){
    if(size=="5mm"){
      EXPERIMENTE<-c("N PROF_PLASMA_NOESY","N PROF_PLASMA_JRES","N PROF_PLASMA_CPMG","N PROF_PLASMA_DIFF","N PROF_PLASMA_DAS_A","N PROF_PLASMA_DAS_L")
    }
    else if(size=="3mm") {
      EXPERIMENTE<-c("N PROF_PLASMA_NOESY_3mm","N PROF_PLASMA_JRES_3mm","N PROF_PLASMA_CPMG_3mm","N PROF_PLASMA_DIFF_3mm","N PROF_PLASMA_DAS_A","N PROF_PLASMA_DAS_L")
    }
  }
  if(solvent=="CSF"){
    if(size=="5mm"){
      EXPERIMENTE<-c("N PROF_PLASMA_NOESY","N PROF_PLASMA_JRES")
    }
    else if(size=="3mm") {
      EXPERIMENTE<-c("N PROF_PLASMA_NOESY_3mm","N PROF_PLASMA_JRES_3mm")
    }
  }
  if(solvent=="MEOH"){
    if(size=="5mm"){
      EXPERIMENTE<-c("N PROF_MEOH_ZG30","N PROF_MEOH_ZGPS","N PROF_PLASMA_NOESY","N PROF_PLASMA_JRES")
    }
    else if(size=="3mm") {
      EXPERIMENTE<-c("N PROF_MEOH_ZG30_3mm","N PROF_MEOH_ZGPS_3mm","N PROF_PLASMA_NOESY_3mm","N PROF_PLASMA_JRES_3mm")
    }
  }
  EXPERIMENTE<-EXPERIMENTE
}

#SCRIPT Plotter----
plotter_v2<-function(data,group,columns,stat){
  data_original<-data
  columns<-c(group,columns)
  data<-select(data,all_of(columns))
  data[data==0] <- NA
  metabolites<-names(data)
  varx<-metabolites[1]
  
  if(stat=="none"){
    plot <- lapply(seq_len(ncol(data)), FUN = function(x) {
      ggplot(data,aes(x = .data[[varx]], y = data[ , x],inherit.aes=FALSE,fill=.data[[varx]],group=.data[[varx]])) +
        geom_boxplot()+
        theme_classic()+
        theme(legend.position="none")+
        stat_boxplot(geom="errorbar",width=0.2)+
        labs(x="", y="", title = colnames(data)[x])
    })
  }
  if(stat=="ttest"){
    plot <- lapply(seq_len(ncol(data)), FUN = function(x) {
      ggplot(data,aes(x = .data[[varx]], y = data[ , x],inherit.aes=FALSE,fill=.data[[varx]],group=.data[[varx]])) +
        geom_boxplot()+
        theme_classic()+
        theme(legend.position="none")+
        labs(x="", y="", title = colnames(data)[x])+
        stat_boxplot(geom="errorbar",width=0.2)+
        stat_compare_means(method="t.test",aes(label=..p.signif..),hide.ns = TRUE)
    })
    
  }
  if(stat=="anova"){
    plot <- lapply(seq_len(ncol(data)), FUN = function(x) {
      ggplot(data,aes(x = .data[[varx]], y = data[ , x],inherit.aes=FALSE,fill=.data[[varx]],group=.data[[varx]])) +
        geom_boxplot()+
        theme_classic()+
        theme(legend.position="none")+
        labs(x="", y="", title = colnames(data)[x])+
        stat_boxplot(geom="errorbar",width=0.2)+
        stat_compare_means(method="anova",aes(label=..p.signif..),hide.ns = TRUE)
    })
  }
  if(stat=="wilcox"){
    plot <- lapply(seq_len(ncol(data)), FUN = function(x) {
      ggplot(data,aes(x = .data[[varx]], y = data[ , x],inherit.aes=FALSE,fill=.data[[varx]],group=.data[[varx]])) +
        geom_boxplot()+
        theme_classic()+
        theme(legend.position="none")+
        labs(x="", y="", title = colnames(data)[x])+
        stat_boxplot(geom="errorbar",width=0.2)+
        stat_compare_means(method="wilcox.text",aes(label=..p.signif..),hide.ns = TRUE)
    })
  }
  if(stat=="kruskal"){
    plot <- lapply(seq_len(ncol(data)), FUN = function(x) {
      ggplot(data,aes(x = .data[[varx]], y = data[ , x],inherit.aes=FALSE,fill=.data[[varx]],group=.data[[varx]])) +
        geom_boxplot()+
        theme_classic()+
        theme(legend.position="none")+
        labs(x="", y="", title = colnames(data)[x])+
        stat_boxplot(geom="errorbar",width=0.2)+
        stat_compare_means(method="kruskal.text",aes(label=..p.signif..),hide.ns = TRUE)
    })
  }
  
  metabolites<-metabolites[-1]
  plot<-plot[-1]
  plot<-plot
  
}
#SCRIPT Normalizer----
normalizor<- function(data,log,center){

  data[data==0] <- NA
  df1<-reshape::melt.data.frame(data)

  if(log=="Yes"){
    df1$value<-log10(df1$value)
  }
  if(center=="Yes"){
    df1$value<-scale(df1$value,center=TRUE, scale=FALSE)
  }

  data2<-pivot_wider(df1,names_from="variable",values_from="value")
}

#SCRIPT Excel Manual----
excel_manual<-function(data,solvent,size,rack,slot,experiments,path){

  #set up the parameters
  SOLVENT<-solvent
  standard_path<-path
  EXPERIMENTE<-experiments

  #create Experiments from input
  if(solvent=="Urine"){
    if(size=="5mm"){
      EXPERIMENTE<-experiments
    }
    else if(size=="3mm"){
      EXPERIMENTE<-experiments
      SOLVENT<-c("Urine_3mm")
    }
  }
  if(solvent=="Urine_Neo"){
    if(size=="5mm"){
      EXPERIMENTE<-experiments
      SOLVENT<-c("Urine")
    }
    else if(size=="3mm"){
      EXPERIMENTE<-experiments
      SOLVENT<-c("Urine_3mm")
    }
  }
  if(solvent=="Plasma"){
    if(size=="5mm"){
      EXPERIMENTE<-experiments}
    else if(size=="3mm") {
      EXPERIMENTE<-experiments
      SOLVENT<-c("Plasma_3mm")
    }
  }
  if(solvent=="CSF"){
    if(size=="5mm"){
      EXPERIMENTE<-experiments
      SOLVENT<-c("CSF")}
    else if(size=="3mm") {
      EXPERIMENTE<-experiments
      SOLVENT<-c("CSF_3mm")
    }
  }
  if(solvent=="MEOH"){
    if(size=="5mm"){
      EXPERIMENTE<-experiments
      SOLVENT<-c("MEOH")}
    else if(size=="3mm") {
      EXPERIMENTE<-experiments
      SOLVENT<-c("MEOH_3mm")
    }
  }
  #read in rest from data
  Date<-Sys.Date()

  #FilePath
  YEAR<-format(as.Date(Date, format="%Y/%m/%d"),"%Y")
  DAY<-format(as.Date(Date, format="%Y/%m/%d"),"%d")
  datapath<-c("data")
  nmr<-c("nmr")
  DISK<-file.path(standard_path,YEAR,datapath,Date,nmr)

  #Names
  NAME<-select(data,Name)
  NAME<- NAME %>%
    rename(NAME = Name)
  TITLE<-NAME
  TITLE<- TITLE %>%
    rename(TITLE = NAME)


  POSI<-slot
  HOLDER<-paste(rack, sprintf("%02d",POSI), sep="")

  #Combine Experiments
  data<-data[order(data$Name),]
  EXPERIMENT<-EXPERIMENTE
  data_final<-cbind(DISK,NAME,SOLVENT,HOLDER,TITLE)
  data_final<-merge(x=data_final,y=EXPERIMENTE)
  data_final<-data_final[order(data_final$HOLDER),]
  data_final<- data_final %>%
    rename(EXPERIMENT = y)
  data_final<-data_final[,c(1,2,3,6,4,5)]
}
#SCRIPT Excel Upload-------------------------------------------
excel_creator<-function(data,solvent,size,rack,experiments,path){

  #set up the parameters
  SOLVENT<-solvent
  standard_path<-path
  EXPERIMENTE<-experiments
  #create Experiments from input
  if(solvent=="Urine"){
    if(size=="5mm"){
      EXPERIMENTE<-experiments
    }
    else if(size=="3mm"){
      EXPERIMENTE<-experiments
      SOLVENT<-c("Urine_3mm")
    }
  }
  if(solvent=="Urine_Neo"){
    if(size=="5mm"){
      EXPERIMENTE<-experiments
      SOLVENT<-c("Urine")
    }
    else if(size=="3mm"){
      EXPERIMENTE<-experiments
      SOLVENT<-c("Urine_3mm")
    }
  }
  if(solvent=="Plasma"){
    if(size=="5mm"){
      EXPERIMENTE<-experiments}
    else if(size=="3mm") {
      EXPERIMENTE<-experiments
      SOLVENT<-c("Plasma_3mm")
    }
  }
  if(solvent=="Media"){
    if(size=="5mm"){
      EXPERIMENTE<-experiments
      SOLVENT<-c("Plasma")}
    else if(size=="3mm") {
      EXPERIMENTE<-experiments
      SOLVENT<-c("Plasma_3mm")
    }
  }
  if(solvent=="CSF"){
    if(size=="5mm"){
      EXPERIMENTE<-experiments
      SOLVENT<-c("CSF")}
    else if(size=="3mm") {
      EXPERIMENTE<-experiments
      SOLVENT<-c("CSF_3mm")
    }
  }
  if(solvent=="MEOH"){
    if(size=="5mm"){
      EXPERIMENTE<-experiments
      SOLVENT<-c("MEOH")}
    else if(size=="3mm") {
      EXPERIMENTE<-experiments
      SOLVENT<-c("MEOH_3mm")
    }
  }
  #read in rest from data
  Date<-Sys.Date()

  #FilePath
  YEAR<-format(as.Date(Date, format="%Y/%m/%d"),"%Y")
  DAY<-format(as.Date(Date, format="%Y/%m/%d"),"%d")
  nmr<-c("nmr")
  datapath<-c("data")
  DISK<-file.path(standard_path,YEAR,datapath,Date,nmr)

  #Names
  NAME<-select(data,Name)
  NAME<- NAME %>%
    rename(NAME = Name)
  TITLE<-NAME
  TITLE<- TITLE %>%
    rename(TITLE = NAME)

  #for sorting purposes only
  order<-nrow(data)
  dataorder<-cbind(data,order)

  #Positions, create 2nd rack if necessary
  if (nrow(data)>192){
    rack_temp2<-data %>%
      slice(193:n())
    data<-data %>%
      slice(1:192)
    rack_temp1<-data %>%
      slice(97:n())
    data<-data %>%
      slice(1:96)
    POSI<-seq(1, nrow(data), by = 1)
    POSI2<-seq(1, nrow(rack_temp1), by = 1)
    POSI3<-seq(1, nrow(rack_temp2), by = 1)
    HOLDER<-paste(rack, sprintf("%02d",POSI), sep="")
    HOLDER2<-paste(rack+1, sprintf("%02d",POSI2), sep="")
    HOLDER3<-paste(rack+2, sprintf("%02d",POSI3), sep="")
    if(rack==4) {
      HOLDER2<-paste(rack+1, sprintf("%02d",POSI2), sep="")
      HOLDER3<-paste(rack-3, sprintf("%02d",POSI3), sep="")
    }
    if(rack==5){
      HOLDER2<-paste(rack-4, sprintf("%02d",POSI2), sep="")
      HOLDER3<-paste(rack-3, sprintf("%02d",POSI3), sep="")
    }
    POSI<-c(POSI,POSI2,POSI3)
    HOLDER<-c(HOLDER,HOLDER2,HOLDER3)

  } else
    if (nrow(data)>96){
      rack_temp<- data %>%
        slice(97:n())
      data<-data %>%
        slice(1:96)
      POSI<-seq(1, nrow(data), by = 1)
      HOLDER<-paste(rack, sprintf("%02d",POSI), sep="")
      POSI2<-seq(1, nrow(rack_temp), by = 1)
      if(rack==5){
        HOLDER2<-paste(1, sprintf("%02d",POSI2), sep="")
      } else {
        HOLDER2<-paste(rack+1, sprintf("%02d",POSI2), sep="")
      }
      POSI<-c(POSI,POSI2)
      HOLDER<-c(HOLDER,HOLDER2)
      #cat(paste("","More than 96 samples (1 rack) were submitted.","A 2nd rack with a higher number (unless 5 then 1) was used from sample 97+",sep="\n"))
    } else {
      POSI<-seq(1, nrow(data), by = 1)
      HOLDER<-paste(rack, sprintf("%02d",POSI), sep="")
    }

  #Combine Experiments
  data<-data[order(data$Name),]
  EXPERIMENT<-EXPERIMENTE
  data_final<-cbind(DISK,NAME,SOLVENT,HOLDER,TITLE)
  data_final<-merge(x=data_final,y=EXPERIMENTE)
  data_final<-data_final[order(data_final$HOLDER),]
  data_final<- data_final %>%
    rename(EXPERIMENT = y)
  # data_final<-sort(data_final$dataorder)
  # data_final=subset(data_final,select=-c(dataorder))
  data_final<-data_final[,c(1,2,3,6,4,5)]
  data_final<-data_final
}
#SCRIPT Extractor----
clean_csv<-function(names,data){

  if(names=="Yes"){
    clean_csv<-clean_names(data)
    if("instrument" %in% colnames(clean_csv)){
      clean_csv<-subset(clean_csv, select = -c(directory,exp_no,proc_no,experiment,pulse_program,aunm,aunmp,instrument,probehead,ns,ds,p1,pld_b_1,pld_b9,rg,swh,td,si,te,phc0,phc1,o1,sf,sr,date,date_2,name,type))
    }
    if("is_b_i_quant_ps" %in% colnames(clean_csv)){
      ID<-subset(clean_csv, select=sample_id)
      ID<-sub("_e.*","",ID [,1])
      clean_csv<-subset(clean_csv, select = -c(measurement_date,reporting_date,is_b_i_quant_ps,is_quant_ps_2_0_0))
      clean_csv<-select(clean_csv, contains("_conc_mmol_l"))
      clean_csv<-clean_csv[seq(1, ncol(clean_csv),3)]
      clean_csv<-cbind(ID,clean_csv)
    }
    if("is_b_i_quant_ur_ne" %in% colnames(clean_csv)){
      clean_csv<-subset(clean_csv, select = -c(measurement_date,reporting_date,is_b_i_quant_ur_ne,is_quant_ur_ne_1_1_0))
      clean_csv<-clean_names(data)
      ID<-subset(clean_csv, select=sample_id)
      ID<-sub("_e.*","",ID [,1])
      crea<-select(clean_csv, contains("creatinine_conc_"))
      clean_csv<-select(clean_csv, contains("_conc_mmol_mol_crea"))
      clean_csv<-cbind(ID,crea,clean_csv)
    }
    if("is_b_i_quant_ur_e" %in% colnames(clean_csv)){
      clean_csv<-subset(clean_csv, select = -c(measurement_date,reporting_date,is_b_i_quant_ur_e,is_quant_ur_e_1_1_0))
      clean_csv<-clean_names(data)
      ID<-subset(clean_csv, select=sample_id)
      ID<-sub("_e.*","",ID [,1])
      crea<-select(clean_csv, contains("creatinine_conc_"))
      clean_csv<-select(clean_csv, contains("_conc_mmol_mol_crea"))
      clean_csv<-cbind(ID,crea,clean_csv)
    }

    #Main Parameters
    names(clean_csv)[names(clean_csv) == 'tptg_mg_d_l'] <- 'TG_mg_dl'
    names(clean_csv)[names(clean_csv) == 'tpch_mg_d_l'] <- 'CHOL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'ldch_mg_d_l'] <- 'LDL-CHOL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'hdch_mg_d_l'] <- 'HDL-CHOL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'tpa1_mg_d_l'] <- 'Apo-A1_mg_dl'
    names(clean_csv)[names(clean_csv) == 'tpa2_mg_d_l'] <- 'Apo-A2_mg_dl'
    names(clean_csv)[names(clean_csv) == 'tpab_mg_d_l'] <- 'Apo-B100_mg_dl'

    #Calculated Figues
    names(clean_csv)[names(clean_csv) == 'ldhd'] <- 'LDL-CHOL_HDL-CHOL'
    names(clean_csv)[names(clean_csv) == 'aba1'] <- 'Apo-B100_Apo-A1'

    #Total Concentration of ApoB carryin particles
    names(clean_csv)[names(clean_csv) == 'tbpn_nmol_l'] <- 'Total_Particles_ApoB_nmol_l'

    #Lipoprotein Main fractions
    names(clean_csv)[names(clean_csv) == 'vlpn_nmol_l'] <- 'VLDL-Particles_nmol_l'
    names(clean_csv)[names(clean_csv) == 'idpn_nmol_l'] <- 'IDL-Particles_nmol_l'
    names(clean_csv)[names(clean_csv) == 'ldpn_nmol_l'] <- 'LDL-Particles_nmol_l'

    #LDL Subfractions
    names(clean_csv)[names(clean_csv) == 'l1pn_nmol_l'] <- 'LDL-1-Particles_nmol_l'
    names(clean_csv)[names(clean_csv) == 'l2pn_nmol_l'] <- 'LDL-2-Particles_nmol_l'
    names(clean_csv)[names(clean_csv) == 'l3pn_nmol_l'] <- 'LDL-3-Particles_nmol_l'
    names(clean_csv)[names(clean_csv) == 'l4pn_nmol_l'] <- 'LDL-4-Particles_nmol_l'
    names(clean_csv)[names(clean_csv) == 'l5pn_nmol_l'] <- 'LDL-5-Particles_nmol_l'
    names(clean_csv)[names(clean_csv) == 'l6pn_nmol_l'] <- 'LDL-6-Particles_nmol_l'

    #Lipoprotein Main fractions
    #TG
    names(clean_csv)[names(clean_csv) == 'vltg_mg_d_l'] <- 'TG_VLDL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'idtg_mg_d_l'] <- 'TG_IDL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'ldtg_mg_d_l'] <- 'TG_LDL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'hdtg_mg_d_l'] <- 'TG_HDL_mg_dl'
    #CHOL
    names(clean_csv)[names(clean_csv) == 'vlch_mg_d_l'] <- 'CHOL_VLDL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'idch_mg_d_l'] <- 'CHOL_IDL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'ldch_mg_d_l'] <- 'CHOL_LDL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'hdch_mg_d_l'] <- 'CHOL_HDL_mg_dl'
    #free CHOL
    names(clean_csv)[names(clean_csv) == 'vlfc_mg_d_l'] <- 'fCHOL_VLDL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'idfc_mg_d_l'] <- 'fCHOL_IDL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'ldfc_mg_d_l'] <- 'fCHOL_LDL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'hdfc_mg_d_l'] <- 'fCHOL_HDL_mg_dl'
    #Phospholipid
    names(clean_csv)[names(clean_csv) == 'vlpl_mg_d_l'] <- 'PHOSL_VLDL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'idpl_mg_d_l'] <- 'PHOSL_IDL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'ldpl_mg_d_l'] <- 'PHOSL_LDL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'hdpl_mg_d_l'] <- 'PHOSL_HDL_mg_dl'

    names(clean_csv)[names(clean_csv) == 'hda1_mg_d_l'] <- 'ApoA1_HDL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'hda2_mg_d_l'] <- 'ApoA2_HDL_mg_dl'

    names(clean_csv)[names(clean_csv) == 'vlab_mg_d_l'] <- 'ApoB_VLDL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'idab_mg_d_l'] <- 'ApoB_IDL_mg_dl'
    names(clean_csv)[names(clean_csv) == 'ldab_mg_d_l'] <- 'ApoB_LDL_mg_dl'

    #VLDL Subfractions
    #TG
    names(clean_csv)[names(clean_csv) == 'v1tg_mg_d_l'] <- 'TG_VLDL-1_mg_dl'
    names(clean_csv)[names(clean_csv) == 'v2tg_mg_d_l'] <- 'TG_VLDL-2_mg_dl'
    names(clean_csv)[names(clean_csv) == 'v3tg_mg_d_l'] <- 'TG_VLDL-3_mg_dl'
    names(clean_csv)[names(clean_csv) == 'v4tg_mg_d_l'] <- 'TG_VLDL-4_mg_dl'
    names(clean_csv)[names(clean_csv) == 'v5tg_mg_d_l'] <- 'TG_VLDL-5_mg_dl'

    #CHOL
    names(clean_csv)[names(clean_csv) == 'v1ch_mg_d_l'] <- 'CHOL_VLDL-1_mg_dl'
    names(clean_csv)[names(clean_csv) == 'v2ch_mg_d_l'] <- 'CHOL_VLDL-2_mg_dl'
    names(clean_csv)[names(clean_csv) == 'v3ch_mg_d_l'] <- 'CHOL_VLDL-3_mg_dl'
    names(clean_csv)[names(clean_csv) == 'v4ch_mg_d_l'] <- 'CHOL_VLDL-4_mg_dl'
    names(clean_csv)[names(clean_csv) == 'v5ch_mg_d_l'] <- 'CHOL_VLDL-5_mg_dl'

    #free CHOL
    names(clean_csv)[names(clean_csv) == 'v1fc_mg_d_l'] <- 'fCHOL_VLDL-1_mg_dl'
    names(clean_csv)[names(clean_csv) == 'v2fc_mg_d_l'] <- 'fCHOL_VLDL-2_mg_dl'
    names(clean_csv)[names(clean_csv) == 'v3fc_mg_d_l'] <- 'fCHOL_VLDL-3_mg_dl'
    names(clean_csv)[names(clean_csv) == 'v4fc_mg_d_l'] <- 'fCHOL_VLDL-4_mg_dl'
    names(clean_csv)[names(clean_csv) == 'v5fc_mg_d_l'] <- 'fCHOL_VLDL-5_mg_dl'

    #PHOSpholipids
    names(clean_csv)[names(clean_csv) == 'v1pl_mg_d_l'] <- 'PHOSL_VLDL-1_mg_dl'
    names(clean_csv)[names(clean_csv) == 'v2pl_mg_d_l'] <- 'PHOSL_VLDL-2_mg_dl'
    names(clean_csv)[names(clean_csv) == 'v3pl_mg_d_l'] <- 'PHOSL_VLDL-3_mg_dl'
    names(clean_csv)[names(clean_csv) == 'v4pl_mg_d_l'] <- 'PHOSL_VLDL-4_mg_dl'
    names(clean_csv)[names(clean_csv) == 'v5pl_mg_d_l'] <- 'PHOSL_VLDL-5_mg_dl'

    #LDL Subfractions
    #TG
    names(clean_csv)[names(clean_csv) == 'l1tg_mg_d_l'] <- 'TG_LDL-1_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l2tg_mg_d_l'] <- 'TG_LDL-2_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l3tg_mg_d_l'] <- 'TG_LDL-3_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l4tg_mg_d_l'] <- 'TG_LDL-4_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l5tg_mg_d_l'] <- 'TG_LDL-5_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l6tg_mg_d_l'] <- 'TG_LDL-6_mg_dl'

    #CHOL
    names(clean_csv)[names(clean_csv) == 'l1ch_mg_d_l'] <- 'CHOL_LDL-1_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l2ch_mg_d_l'] <- 'CHOL_LDL-2_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l3ch_mg_d_l'] <- 'CHOL_LDL-3_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l4ch_mg_d_l'] <- 'CHOL_LDL-4_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l5ch_mg_d_l'] <- 'CHOL_LDL-5_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l6ch_mg_d_l'] <- 'CHOL_LDL-6_mg_dl'

    #free CHOL
    names(clean_csv)[names(clean_csv) == 'l1fc_mg_d_l'] <- 'fCHOL_LDL-1_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l2fc_mg_d_l'] <- 'fCHOL_LDL-2_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l3fc_mg_d_l'] <- 'fCHOL_LDL-3_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l4fc_mg_d_l'] <- 'fCHOL_LDL-4_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l5fc_mg_d_l'] <- 'fCHOL_LDL-5_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l6fc_mg_d_l'] <- 'fCHOL_LDL-6_mg_dl'

    #PHOSpholipids
    names(clean_csv)[names(clean_csv) == 'l1pl_mg_d_l'] <- 'PHOSL_LDL-1_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l2pl_mg_d_l'] <- 'PHOSL_LDL-2_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l3pl_mg_d_l'] <- 'PHOSL_LDL-3_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l4pl_mg_d_l'] <- 'PHOSL_LDL-4_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l5pl_mg_d_l'] <- 'PHOSL_LDL-5_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l6pl_mg_d_l'] <- 'PHOSL_LDL-6_mg_dl'

    #ApoB
    names(clean_csv)[names(clean_csv) == 'l1ab_mg_d_l'] <- 'ApoB_LDL-1_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l2ab_mg_d_l'] <- 'ApoB_LDL-2_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l3ab_mg_d_l'] <- 'ApoB_LDL-3_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l4ab_mg_d_l'] <- 'ApoB_LDL-4_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l5ab_mg_d_l'] <- 'ApoB_LDL-5_mg_dl'
    names(clean_csv)[names(clean_csv) == 'l6ab_mg_d_l'] <- 'ApoB_LDL-6_mg_dl'

    #HDL Subfractions
    #TG
    names(clean_csv)[names(clean_csv) == 'h1tg_mg_d_l'] <- 'TG_HDL-1_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h2tg_mg_d_l'] <- 'TG_HDL-2_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h3tg_mg_d_l'] <- 'TG_HDL-3_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h4tg_mg_d_l'] <- 'TG_HDL-4_mg_dl'

    #CHOL
    names(clean_csv)[names(clean_csv) == 'h1ch_mg_d_l'] <- 'CHOL_HDL-1_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h2ch_mg_d_l'] <- 'CHOL_HDL-2_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h3ch_mg_d_l'] <- 'CHOL_HDL-3_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h4ch_mg_d_l'] <- 'CHOL_HDL-4_mg_dl'

    #free CHOL
    names(clean_csv)[names(clean_csv) == 'h1fc_mg_d_l'] <- 'fCHOL_HDL-1_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h2fc_mg_d_l'] <- 'fCHOL_HDL-2_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h3fc_mg_d_l'] <- 'fCHOL_HDL-3_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h4fc_mg_d_l'] <- 'fCHOL_HDL-4_mg_dl'

    #Phospholipids
    names(clean_csv)[names(clean_csv) == 'h1pl_mg_d_l'] <- 'PHOSL_HDL-1_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h2pl_mg_d_l'] <- 'PHOSL_HDL-2_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h3pl_mg_d_l'] <- 'PHOSL_HDL-3_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h4pl_mg_d_l'] <- 'PHOSL_HDL-4_mg_dl'

    #ApoA1
    names(clean_csv)[names(clean_csv) == 'h1a1_mg_d_l'] <- 'ApoA1_HDL-1_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h2a1_mg_d_l'] <- 'ApoA1_HDL-2_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h3a1_mg_d_l'] <- 'ApoA1_HDL-3_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h4a1_mg_d_l'] <- 'ApoA1_HDL-4_mg_dl'

    #ApoA2
    names(clean_csv)[names(clean_csv) == 'h1a2_mg_d_l'] <- 'ApoA2_HDL-1_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h2a2_mg_d_l'] <- 'ApoA2_HDL-2_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h3a2_mg_d_l'] <- 'ApoA2_HDL-3_mg_dl'
    names(clean_csv)[names(clean_csv) == 'h4a2_mg_d_l'] <- 'ApoA2_HDL-4_mg_dl'
    #end----

    #plasma names----
    names(clean_csv)[names(clean_csv) == 'ethanol_conc_mmol_l'] <- 'Ethanol'
    names(clean_csv)[names(clean_csv) == 'trimethylamine_n_oxide_conc_mmol_l'] <- 'Trimethylamine N-oxide'
    names(clean_csv)[names(clean_csv) == 'x2_aminobutyric_acid_conc_mmol_l'] <- '2-Aminobutyrate'
    names(clean_csv)[names(clean_csv) == 'alanine_conc_mmol_l'] <- 'Alanine'
    names(clean_csv)[names(clean_csv) == 'asparagine_conc_mmol_l'] <- 'Asparagine'
    names(clean_csv)[names(clean_csv) == 'creatine_conc_mmol_l'] <- 'Creatine'
    names(clean_csv)[names(clean_csv) == 'creatinine_conc_mmol_l'] <- 'Creatinine'
    names(clean_csv)[names(clean_csv) == 'glutamic_acid_conc_mmol_l'] <- 'Glutamate'
    names(clean_csv)[names(clean_csv) == 'glutamine_conc_mmol_l'] <- 'Glutamine'
    names(clean_csv)[names(clean_csv) == 'glycine_conc_mmol_l'] <- 'Glycine'
    names(clean_csv)[names(clean_csv) == 'histidine_conc_mmol_l'] <- 'Histidine'
    names(clean_csv)[names(clean_csv) == 'isoleucine_conc_mmol_l'] <- 'Isoleucine'
    names(clean_csv)[names(clean_csv) == 'leucine_conc_mmol_l'] <- 'Leucine'
    names(clean_csv)[names(clean_csv) == 'lysine_conc_mmol_l'] <- 'Lysine'
    names(clean_csv)[names(clean_csv) == 'methionine_conc_mmol_l'] <- 'Methionine'
    names(clean_csv)[names(clean_csv) == 'n_n_dimethylglycine_conc_mmol_l'] <- 'N,N-Dimethylglycine'
    names(clean_csv)[names(clean_csv) == 'ornithine_conc_mmol_l'] <- 'Ornithine'
    names(clean_csv)[names(clean_csv) == 'phenylalanine_conc_mmol_l'] <- 'Phenylalanine'
    names(clean_csv)[names(clean_csv) == 'proline_conc_mmol_l'] <- 'Proline'
    names(clean_csv)[names(clean_csv) == 'sarcosine_conc_mmol_l'] <- 'Sarcosine'
    names(clean_csv)[names(clean_csv) == 'threonine_conc_mmol_l'] <- 'Threonine'
    names(clean_csv)[names(clean_csv) == 'tyrosine_conc_mmol_l'] <- 'Tyrosine'
    names(clean_csv)[names(clean_csv) == 'valine_conc_mmol_l'] <- 'Valine'
    names(clean_csv)[names(clean_csv) == 'x2_hydroxybutyric_acid_conc_mmol_l'] <- '2-Hydroxybutyrate'
    names(clean_csv)[names(clean_csv) == 'acetic_acid_conc_mmol_l'] <- 'Acetate'
    names(clean_csv)[names(clean_csv) == 'citric_acid_conc_mmol_l'] <- 'Citrate'
    names(clean_csv)[names(clean_csv) == 'formic_acid_conc_mmol_l'] <- 'Formate'
    names(clean_csv)[names(clean_csv) == 'lactic_acid_conc_mmol_l'] <- 'Lactate'
    names(clean_csv)[names(clean_csv) == 'succinic_acid_conc_mmol_l'] <- 'Succinate'
    names(clean_csv)[names(clean_csv) == 'choline_conc_mmol_l'] <- 'Choline'
    names(clean_csv)[names(clean_csv) == 'x2_oxoglutaric_acid_conc_mmol_l'] <- '2-Oxoglutarate'
    names(clean_csv)[names(clean_csv) == 'x3_hydroxybutyric_acid_conc_mmol_l'] <- '3-Hydroxybutyrate'
    names(clean_csv)[names(clean_csv) == 'acetoacetic_acid_conc_mmol_l'] <- 'Acetoacetate'
    names(clean_csv)[names(clean_csv) == 'acetone_conc_mmol_l'] <- 'Acetone'
    names(clean_csv)[names(clean_csv) == 'pyruvic_acid_conc_mmol_l'] <- 'Pyruvate'
    names(clean_csv)[names(clean_csv) == 'd_galactose_conc_mmol_l'] <- 'Galactose'
    names(clean_csv)[names(clean_csv) == 'glucose_conc_mmol_l'] <- 'Glucose'
    names(clean_csv)[names(clean_csv) == 'glycerol_conc_mmol_l'] <- 'Glycerol'
    names(clean_csv)[names(clean_csv) == 'dimethylsulfone_conc_mmol_l'] <- 'Dimethylsufone'
    names(clean_csv)[names(clean_csv) == 'ca_edta_conc_mmol_l'] <- 'Ca_EDTA'
    names(clean_csv)[names(clean_csv) == 'k_edta_conc_mmol_l'] <- 'K_EDTA'

    #urine names ----
    names(clean_csv)[names(clean_csv) == 'creatinine_conc_mmol_l'] <- 'Creatinine'
    names(clean_csv)[names(clean_csv) == 'ethanol_conc_mmol_mol_crea'] <- 'Ethanol'
    names(clean_csv)[names(clean_csv) == 'isopropanol_conc_mmol_mol_crea'] <- 'Isopropanol'
    names(clean_csv)[names(clean_csv) == 'methanol_conc_mmol_mol_crea'] <- 'Methanol'
    names(clean_csv)[names(clean_csv) == 'propylene_glycol_conc_mmol_mol_crea'] <- 'Propylene glycol'
    names(clean_csv)[names(clean_csv) == 'x1_methylguanidine_conc_mmol_mol_crea'] <- 'Methylguanidine'
    names(clean_csv)[names(clean_csv) == 'dimethylamine_conc_mmol_mol_crea'] <- 'Dimethylamine'
    names(clean_csv)[names(clean_csv) == 'trimethylamine_conc_mmol_mol_crea'] <- 'Trimethylamine'
    names(clean_csv)[names(clean_csv) == 'tyramine_conc_mmol_mol_crea'] <- 'Tyramine'
    names(clean_csv)[names(clean_csv) == 'x1_methylhistidine_conc_mmol_mol_crea'] <- '1-Methylhistidine'
    names(clean_csv)[names(clean_csv) == 'x2_furoylglycine_conc_mmol_mol_crea'] <- '2-Furoylglycine'
    names(clean_csv)[names(clean_csv) == 'x3_aminoisobutyric_acid_conc_mmol_mol_crea'] <- '3-Aminoisobutyrate'
    names(clean_csv)[names(clean_csv) == 'x3_methylcrotonylglycine_conc_mmol_mol_crea'] <- '3-Methylcrotonylglycine'
    names(clean_csv)[names(clean_csv) == 'x4_aminobutyric_acid_conc_mmol_mol_crea'] <- '4-Aminobutyrate'
    names(clean_csv)[names(clean_csv) == 'x5_aminopentanoic_acid_conc_mmol_mol_crea'] <- '5-Aminopentanoate'
    names(clean_csv)[names(clean_csv) == 'alanine_conc_mmol_mol_crea'] <- 'Alanine'
    names(clean_csv)[names(clean_csv) == 'arginine_conc_mmol_mol_crea'] <- 'Arginine'
    names(clean_csv)[names(clean_csv) == 'argininosuccinic_acid_conc_mmol_mol_crea'] <- 'Argininosuccinate'
    names(clean_csv)[names(clean_csv) == 'betaine_conc_mmol_mol_crea'] <- 'Betaine'
    names(clean_csv)[names(clean_csv) == 'citrulline_conc_mmol_mol_crea'] <- 'Citrulline'
    names(clean_csv)[names(clean_csv) == 'creatine_conc_mmol_mol_crea'] <- 'Creatine'
    names(clean_csv)[names(clean_csv) == 'cystine_conc_mmol_mol_crea'] <- 'Cystine'
    names(clean_csv)[names(clean_csv) == 'dl_alloisoleucine_conc_mmol_mol_crea'] <- 'DL-Alloisoleucine'
    names(clean_csv)[names(clean_csv) == 'dl_tyrosine_conc_mmol_mol_crea'] <- 'Tyrosine'
    names(clean_csv)[names(clean_csv) == 'glutamic_acid_conc_mmol_mol_crea'] <- 'Glutamate'
    names(clean_csv)[names(clean_csv) == 'glutamine_conc_mmol_mol_crea'] <- 'Glutamine'
    names(clean_csv)[names(clean_csv) == 'glycine_conc_mmol_mol_crea'] <- 'Glycine'
    names(clean_csv)[names(clean_csv) == 'guanidinoacetic_acid_conc_mmol_mol_crea'] <- 'Guanidinoacetate'
    names(clean_csv)[names(clean_csv) == 'isobutyrylglycine_conc_mmol_mol_crea'] <- 'Isobutyrylglycine'
    names(clean_csv)[names(clean_csv) == 'l_carnosine_conc_mmol_mol_crea'] <- 'L-Carnosine'
    names(clean_csv)[names(clean_csv) == 'l_homocystine_conc_mmol_mol_crea'] <- 'L-Homocystine'
    names(clean_csv)[names(clean_csv) == 'l_isoleucine_conc_mmol_mol_crea'] <- 'L-Isoleucine'
    names(clean_csv)[names(clean_csv) == 'l_pyroglutamic_acid_conc_mmol_mol_crea'] <- 'L-Pyroglutamate'
    names(clean_csv)[names(clean_csv) == 'l_tryptophan_conc_mmol_mol_crea'] <- 'L-Tryptophan'
    names(clean_csv)[names(clean_csv) == 'leucine_conc_mmol_mol_crea'] <- 'Leucine'
    names(clean_csv)[names(clean_csv) == 'methionine_conc_mmol_mol_crea'] <- 'Methionine'
    names(clean_csv)[names(clean_csv) == 'n_n_dimethylglycine_conc_mmol_mol_crea'] <- 'N,N-Dimethylglycine'
    names(clean_csv)[names(clean_csv) == 'n_acetylaspartic_acid_conc_mmol_mol_crea'] <- 'N-Acetylaspartate'
    names(clean_csv)[names(clean_csv) == 'n_acetylglutamate_conc_mmol_mol_crea'] <- 'N-Acetylglutamate'
    names(clean_csv)[names(clean_csv) == 'n_acetylphenylalanine_conc_mmol_mol_crea'] <- 'N-Acetylphenylalanine'
    names(clean_csv)[names(clean_csv) == 'n_acetyltyrosine_conc_mmol_mol_crea'] <- 'N-Acetyltyrosine'
    names(clean_csv)[names(clean_csv) == 'n_isovaleroylglycine_conc_mmol_mol_crea'] <- 'N-Isovaleroylglycine'
    names(clean_csv)[names(clean_csv) == 'phenylalanine_conc_mmol_mol_crea'] <- 'Phenylalanine'
    names(clean_csv)[names(clean_csv) == 'proline_betaine_conc_mmol_mol_crea'] <- 'Proline betaine'
    names(clean_csv)[names(clean_csv) == 'propionylglycine_conc_mmol_mol_crea'] <- 'Propionylglycine'
    names(clean_csv)[names(clean_csv) == 'sarcosine_conc_mmol_mol_crea'] <- 'Sarcosine'
    names(clean_csv)[names(clean_csv) == 'taurine_conc_mmol_mol_crea'] <- 'Taurine'
    names(clean_csv)[names(clean_csv) == 'tiglylglycine_conc_mmol_mol_crea'] <- 'Tiglylglycine'
    names(clean_csv)[names(clean_csv) == 'valine_conc_mmol_mol_crea'] <- 'Valine'
    names(clean_csv)[names(clean_csv) == 'x2_hydroxyphenylacetic_acid_conc_mmol_mol_crea'] <- '2-Hydroxyphenylacetate'
    names(clean_csv)[names(clean_csv) == 'x3_phenyllactic_acid_conc_mmol_mol_crea'] <- '3-Phenyllactate'
    names(clean_csv)[names(clean_csv) == 'x4_aminohippuric_acid_conc_mmol_mol_crea'] <- '4-Aminohippurate'
    names(clean_csv)[names(clean_csv) == 'x4_ethylphenol_conc_mmol_mol_crea'] <- '4-Ethylphenol'
    names(clean_csv)[names(clean_csv) == 'x4_hydroxyhippuric_acid_conc_mmol_mol_crea'] <- '4-Hydroxyhippurate'
    names(clean_csv)[names(clean_csv) == 'x4_hydroxyphenylacetic_acid_conc_mmol_mol_crea'] <- '4-Hydroxyphenylacetate'
    names(clean_csv)[names(clean_csv) == 'x4_hydroxyphenyllactic_acid_conc_mmol_mol_crea'] <- '4-Hyxroxyphenyllactate'
    names(clean_csv)[names(clean_csv) == 'benzoic_acid_conc_mmol_mol_crea'] <- 'Benzoate'
    names(clean_csv)[names(clean_csv) == 'd_mandelic_acid_conc_mmol_mol_crea'] <- 'D-Mandelicate'
    names(clean_csv)[names(clean_csv) == 'hippuric_acid_conc_mmol_mol_crea'] <- 'Hippurate'
    names(clean_csv)[names(clean_csv) == 'phenylacetic_acid_conc_mmol_mol_crea'] <- 'Phenylacetate'
    names(clean_csv)[names(clean_csv) == 'phenylpyruvic_acid_conc_mmol_mol_crea'] <- 'Phenylpyruvate'
    names(clean_csv)[names(clean_csv) == 'pyrocatechol_conc_mmol_mol_crea'] <- 'Pyrocatechol'
    names(clean_csv)[names(clean_csv) == 'syringic_acid_conc_mmol_mol_crea'] <- 'Syringate'
    names(clean_csv)[names(clean_csv) == 'x5_aminolevulinic_acid_conc_mmol_mol_crea'] <- '5-Aminolevulinate'
    names(clean_csv)[names(clean_csv) == 'acetic_acid_conc_mmol_mol_crea'] <- 'Acetate'
    names(clean_csv)[names(clean_csv) == 'citric_acid_conc_mmol_mol_crea'] <- 'Citrate'
    names(clean_csv)[names(clean_csv) == 'e_glutaconic_acid_conc_mmol_mol_crea'] <- 'Glutaconate'
    names(clean_csv)[names(clean_csv) == 'ethylmalonic_acid_conc_mmol_mol_crea'] <- 'Ethylmalonate'
    names(clean_csv)[names(clean_csv) == 'formic_acid_conc_mmol_mol_crea'] <- 'Formate'
    names(clean_csv)[names(clean_csv) == 'fumaric_acid_conc_mmol_mol_crea'] <- 'Fumarate'
    names(clean_csv)[names(clean_csv) == 'glutaric_acid_conc_mmol_mol_crea'] <- 'Glutarate'
    names(clean_csv)[names(clean_csv) == 'imidazole_conc_mmol_mol_crea'] <- 'Imidazole'
    names(clean_csv)[names(clean_csv) == 'lactic_acid_conc_mmol_mol_crea'] <- 'Lactate'
    names(clean_csv)[names(clean_csv) == 'maleic_acid_conc_mmol_mol_crea'] <- 'Maleic acid'
    names(clean_csv)[names(clean_csv) == 'methylmalonic_acid_conc_mmol_mol_crea'] <- 'Methylmalonate'
    names(clean_csv)[names(clean_csv) == 'propionic_acid_conc_mmol_mol_crea'] <- 'Propionate'
    names(clean_csv)[names(clean_csv) == 'succinic_acid_conc_mmol_mol_crea'] <- 'Succinate'
    names(clean_csv)[names(clean_csv) == 'tartaric_acid_conc_mmol_mol_crea'] <- 'Tartarate'
    names(clean_csv)[names(clean_csv) == 'trigonelline_conc_mmol_mol_crea'] <- 'Trigonelline'
    names(clean_csv)[names(clean_csv) == 'xanthurenic_acid_conc_mmol_mol_crea'] <- 'Xanthurenate'
    names(clean_csv)[names(clean_csv) == 'choline_conc_mmol_mol_crea'] <- 'Choline'
    names(clean_csv)[names(clean_csv) == 'd_panthenol_conc_mmol_mol_crea'] <- 'D-Panthenol'
    names(clean_csv)[names(clean_csv) == 'l_ascorbic_acid_conc_mmol_mol_crea'] <- 'L-Ascorbate'
    names(clean_csv)[names(clean_csv) == 'pantothenic_acid_conc_mmol_mol_crea'] <- 'Pantothenate'
    names(clean_csv)[names(clean_csv) == 'paracetamol_conc_mmol_mol_crea'] <- 'Paracetamol'
    names(clean_csv)[names(clean_csv) == 'paracetamol_glucuronide_conc_mmol_mol_crea'] <- 'Paracetamol-Glucuronide'
    names(clean_csv)[names(clean_csv) == 'x2_hydroxy_4_methylvaleric_acid_conc_mmol_mol_crea'] <- '2-Hydroxy-4-Methylvalerate'
    names(clean_csv)[names(clean_csv) == 'x2_hydroxyisovaleric_acid_conc_mmol_mol_crea'] <- '2-Hydroxyisovalerate'
    names(clean_csv)[names(clean_csv) == 'x2_methylsuccinic_acid_conc_mmol_mol_crea'] <- '2-Methylsuccinate'
    names(clean_csv)[names(clean_csv) == 'x3_hydroxy_3_methylglutaric_acid_conc_mmol_mol_crea'] <- '3-Hydroxy-3-Methylglutarate'
    names(clean_csv)[names(clean_csv) == 'x3_hydroxyisovaleric_acid_conc_mmol_mol_crea'] <- '3-Hydroxyisovalerate'
    names(clean_csv)[names(clean_csv) == 'x3_hydroxyvaleric_acid_conc_mmol_mol_crea'] <- '3-Hydroxyvalerate'
    names(clean_csv)[names(clean_csv) == 'x3_methylglutaconic_acid_conc_mmol_mol_crea'] <- '3-Methylglutaconate'
    names(clean_csv)[names(clean_csv) == 'butyric_acid_conc_mmol_mol_crea'] <- 'Butyrate'
    names(clean_csv)[names(clean_csv) == 'citraconic_acid_conc_mmol_mol_crea'] <- 'Citraconate'
    names(clean_csv)[names(clean_csv) == 'l_citramalic_acid_conc_mmol_mol_crea'] <- 'L-Citramalate'
    names(clean_csv)[names(clean_csv) == 'pimelic_acid_conc_mmol_mol_crea'] <- 'Pimelate'
    names(clean_csv)[names(clean_csv) == 'thymol_conc_mmol_mol_crea'] <- 'Thymol'
    names(clean_csv)[names(clean_csv) == 'x3_hydroxyglutaric_acid_conc_mmol_mol_crea'] <- '3-Hydroxyglutarate'
    names(clean_csv)[names(clean_csv) == 'x3_hydroxypropionic_acid_conc_mmol_mol_crea'] <- '3-Hydroxypropionate'
    names(clean_csv)[names(clean_csv) == 'd_galactonic_acid_conc_mmol_mol_crea'] <- 'D-Galactonate'
    names(clean_csv)[names(clean_csv) == 'd_gluconic_acid_conc_mmol_mol_crea'] <- 'D-Gluconate'
    names(clean_csv)[names(clean_csv) == 'glycolic_acid_conc_mmol_mol_crea'] <- 'Glycolate'
    names(clean_csv)[names(clean_csv) == 'malic_acid_conc_mmol_mol_crea'] <- 'Malate'
    names(clean_csv)[names(clean_csv) == 'x2_ketobutyric_acid_conc_mmol_mol_crea'] <- '2-Ketobutyrate'
    names(clean_csv)[names(clean_csv) == 'x2_oxoglutaric_acid_conc_mmol_mol_crea'] <- '2-Oxoglutarate'
    names(clean_csv)[names(clean_csv) == 'x2_oxoisocaproic_acid_conc_mmol_mol_crea'] <- '2-Oxoisocaproate'
    names(clean_csv)[names(clean_csv) == 'x2_oxoisovaleric_acid_conc_mmol_mol_crea'] <- '2-Oxoisovalerate'
    names(clean_csv)[names(clean_csv) == 'x3_hydroxybutyric_acid_conc_mmol_mol_crea'] <- '3-Hydroxybutyrate'
    names(clean_csv)[names(clean_csv) == 'x3_methyl_2_oxovaleric_acid_conc_mmol_mol_crea'] <- '3-Methyl-2-Oxovalerate'
    names(clean_csv)[names(clean_csv) == 'x4_hydroxyphenylpyruvic_acid_conc_mmol_mol_crea'] <- '4-Hydroxyphenylpyruvate'
    names(clean_csv)[names(clean_csv) == 'acetoacetic_acid_conc_mmol_mol_crea'] <- 'Acetoacetate'
    names(clean_csv)[names(clean_csv) == 'acetoine_conc_mmol_mol_crea'] <- 'Acetoine'
    names(clean_csv)[names(clean_csv) == 'acetone_conc_mmol_mol_crea'] <- 'Acetone'
    names(clean_csv)[names(clean_csv) == 'dl_kynurenin_conc_mmol_mol_crea'] <- 'DL-Kynurenin'
    names(clean_csv)[names(clean_csv) == 'oxaloacetic_acid_conc_mmol_mol_crea'] <- 'Oxaloacetate'
    names(clean_csv)[names(clean_csv) == 'pyruvic_acid_conc_mmol_mol_crea'] <- 'Pyruvate'
    names(clean_csv)[names(clean_csv) == 'succinylacetone_conc_mmol_mol_crea'] <- 'Succinylacetone'
    names(clean_csv)[names(clean_csv) == 'x1_3_dimethyluric_acid_conc_mmol_mol_crea'] <- '1,3-Dimethylurate'
    names(clean_csv)[names(clean_csv) == 'x1_methyladenosine_conc_mmol_mol_crea'] <- '1-Methyladenosine'
    names(clean_csv)[names(clean_csv) == 'x1_methylhydantoin_conc_mmol_mol_crea'] <- '1-Methylhydantoin'
    names(clean_csv)[names(clean_csv) == 'x1_methylnicotinamide_conc_mmol_mol_crea'] <- '1-Methylnicotinamie'
    names(clean_csv)[names(clean_csv) == 'x4_pyridoxic_acid_conc_mmol_mol_crea'] <- '4-Pyridoxate'
    names(clean_csv)[names(clean_csv) == 'adenine_conc_mmol_mol_crea'] <- 'Adenine'
    names(clean_csv)[names(clean_csv) == 'adenosine_conc_mmol_mol_crea'] <- 'Adenosine'
    names(clean_csv)[names(clean_csv) == 'allantoin_conc_mmol_mol_crea'] <- 'Allantoin'
    names(clean_csv)[names(clean_csv) == 'allopurinol_conc_mmol_mol_crea'] <- 'Allopurinol'
    names(clean_csv)[names(clean_csv) == 'caffeine_conc_mmol_mol_crea'] <- 'Caffeine'
    names(clean_csv)[names(clean_csv) == 'cytosine_conc_mmol_mol_crea'] <- 'Cytosine'
    names(clean_csv)[names(clean_csv) == 'dihydrothymine_conc_mmol_mol_crea'] <- 'Dihydrothymine'
    names(clean_csv)[names(clean_csv) == 'dihydrouracil_conc_mmol_mol_crea'] <- 'Dihydrouracil'
    names(clean_csv)[names(clean_csv) == 'inosine_conc_mmol_mol_crea'] <- 'Inosine'
    names(clean_csv)[names(clean_csv) == 'neopterin_conc_mmol_mol_crea'] <- 'Neopterin'
    names(clean_csv)[names(clean_csv) == 'orotic_acid_conc_mmol_mol_crea'] <- 'Orotate'
    names(clean_csv)[names(clean_csv) == 'oxypurinol_conc_mmol_mol_crea'] <- 'Oxypurinol'
    names(clean_csv)[names(clean_csv) == 'quinolinic_acid_conc_mmol_mol_crea'] <- 'Quinolinate'
    names(clean_csv)[names(clean_csv) == 'theobromine_conc_mmol_mol_crea'] <- 'Theobromine'
    names(clean_csv)[names(clean_csv) == 'thymine_conc_mmol_mol_crea'] <- 'Thymine'
    names(clean_csv)[names(clean_csv) == 'uracil_conc_mmol_mol_crea'] <- 'Uracil'
    names(clean_csv)[names(clean_csv) == 'uridine_conc_mmol_mol_crea'] <- 'Uridine'
    names(clean_csv)[names(clean_csv) == 'd_galactose_conc_mmol_mol_crea'] <- 'D-Galactose'
    names(clean_csv)[names(clean_csv) == 'd_glucose_conc_mmol_mol_crea'] <- 'D-Glucose'
    names(clean_csv)[names(clean_csv) == 'd_lactose_conc_mmol_mol_crea'] <- 'D-Lactose'
    names(clean_csv)[names(clean_csv) == 'd_mannitol_conc_mmol_mol_crea'] <- 'D-Mannitol'
    names(clean_csv)[names(clean_csv) == 'd_mannose_conc_mmol_mol_crea'] <- 'D-Mannose'
    names(clean_csv)[names(clean_csv) == 'galactitol_conc_mmol_mol_crea'] <- 'Galactitol'
    names(clean_csv)[names(clean_csv) == 'glycerol_conc_mmol_mol_crea'] <- 'Glycerol'
    names(clean_csv)[names(clean_csv) == 'l_fucose_conc_mmol_mol_crea'] <- 'L-Fucose'
    names(clean_csv)[names(clean_csv) == 'l_threonic_acid_conc_mmol_mol_crea'] <- 'L-Threonate'
    names(clean_csv)[names(clean_csv) == 'd_xylose_conc_mmol_mol_crea'] <- 'D-Xylose'
    names(clean_csv)[names(clean_csv) == 'myo_inositol_conc_mmol_mol_crea'] <- 'Myo-Inositol'



    clean_csv<-clean_csv
  }
  #rest----
  # #Only keeps concentration columns
  # clean_csv=data[keeps]
  # clean_csv<-clean_names(clean_csv)
  # #cleans up all spaces and puts _ in there
  # newnames<-sub("_conc.*","",colnames(clean_csv))
  # #removes everything after _conc
  #oldnames<-colnames(clean_csv)
  # #gets the old names into a string
  #clean_csv<- clean_csv %>%
  #   rename_with(~ newnames[which(oldnames == .x)], .cols = oldnames)
  # #replaces oldnames with newnames
  # newnames2<-sub("_e.*","",clean_csv [,1])
  # clean_csv$sample_id<-newnames2


  #replaces the Sample names with just their Code
  #  write.csv(clean_csv, paste0(format(Sys.time(), "%d-%b-%Y %H.%M.%S"), ".csv"),row.names=FALSE)
  #  return("Bruker CSV cleaned and resaved with timestamp.")

  if(names=="No"){
    clean_csv<-clean_names(data)

    if("instrument" %in% colnames(clean_csv)){
      clean_csv<-subset(clean_csv, select = -c(directory,exp_no,proc_no,experiment,pulse_program,aunm,aunmp,instrument,probehead,ns,ds,p1,pld_b_1,pld_b9,rg,swh,td,si,te,phc0,phc1,o1,sf,sr,date,date_2,name,type))
      clean_csv<-clean_csv
    }
    if("is_b_i_quant_ur_ne" %in% colnames(clean_csv)){
      ID<-subset(clean_csv, select=sample_id)
      ID<-sub("_e.*","",ID [,1])
      clean_csv<-subset(clean_csv, select = -c(measurement_date,reporting_date,is_b_i_quant_ur_ne,is_quant_ur_ne_1_1_0))
      crea<-select(clean_csv, contains("creatinine_conc_"))
      clean_csv<-select(clean_csv, contains("_conc_mmol_mol_crea"))
      clean_csv<-clean_csv[seq(1, ncol(clean_csv),3)]
      clean_csv<-cbind(ID,crea,clean_csv)
    }
    if("is_b_i_quant_ur_e" %in% colnames(clean_csv)){
      ID<-subset(clean_csv, select=sample_id)
      ID<-sub("_e.*","",ID [,1])
      clean_csv<-subset(clean_csv, select = -c(measurement_date,reporting_date,is_b_i_quant_ur_e,is_quant_ur_e_1_1_0))
      crea<-select(clean_csv, contains("creatinine_conc_"))
      clean_csv<-select(clean_csv, contains("_conc_mmol_mol_crea"))
      clean_csv<-clean_csv[seq(1, ncol(clean_csv),3)]
      clean_csv<-cbind(ID,crea,clean_csv)
    }
    if("is_b_i_quant_ps" %in% colnames(clean_csv)){
      ID<-subset(clean_csv, select=sample_id)
      ID<-sub("_e.*","",ID [,1])
      clean_csv<-subset(clean_csv, select = -c(measurement_date,reporting_date,is_b_i_quant_ps,is_quant_ps_2_0_0))
      clean_csv<-select(clean_csv, contains("_conc_mmol_l"))
      clean_csv<-clean_csv[seq(1, ncol(clean_csv),3)]
      clean_csv<-cbind(ID,clean_csv)
    }

  }
  clean_csv<-clean_csv
}

  #Waiting screen----
  waiting_screen <- tagList(
                            spin_hexdots(),br(),br(),br(),br(),
                            h4(style="color:black","Generating Plots...")
  )
# UI ----
app_ui <-dashboardPage(
  dashboardHeader(title="MoonNMR"),
  dashboardSidebar(
     sidebarMenu(collapsed=FALSE,
                 menuItem("Home", tabName = "home", icon = icon("home"),selected=T),
                 menuItem("Preparation",icon=icon("dolly"),
                          menuSubItem("IVDr Upload", tabName = "creator", icon = icon("truck")),
                          menuSubItem("IVDr Manual", tabName = "manual", icon = icon("people-carry")),
                          menuSubItem("ICON Upload", tabName = "iconup", icon = icon("car-side")),
                          menuSubItem("ICON Manual", tabName = "iconman", icon = icon("wrench"))
                 ),
                 menuItem("Extraction", icon=icon("bacon"),
                          menuSubItem("XML", tabName = "xml", icon = icon("file-import")),
                          menuSubItem("Bruker File", tabName = "extractor", icon = icon("crow"))
                 ),
                 menuItem("Processing", icon=icon("table"),
                          menuSubItem("Delete",tabName="deleter",icon=icon("trash")),
                          menuSubItem("Combine",tabName="combiner",icon=icon("handshake"))
                 ),
                 menuItem("Analysis",icon=icon("flask"),
                          menuSubItem("Normalize", tabName = "normalizer", icon = icon("balance-scale-right")),
                          menuSubItem("Bar Plots", tabName = "plots", icon = icon("chart-bar")),
                          menuSubItem("PCA", tabName = "pca", icon = icon("chart-area"))
                 )
    )
  ),
  dashboardBody(

    tags$head(tags$style(HTML('
            /* logo */
        .skin-blue .main-header .logo {
                              background-color: #001787;
                              color:#c4ff00;
         /* logo when hovered */
        .skin-blue .main-header .logo:hover {
                              background-color: #8d89ad;
                              }
                              }
        /* navbar (rest of the header) */
        .skin-blue .main-header .navbar {
                              background-color: #001787;
                              color:#c4ff00;
        }
        .skin-blue .main-header .sidebar-toffle {
                              background-color: #001787;
                              color:#c4ff00;
                              }

        /* main sidebar */
        .skin-blue .main-sidebar {
                              background-color: #001787;
                              color:#c4ff00;
        }
        /* active selected tab in the sidebarmenu */
        .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                              background-color: #001787;
                              color:#c4ff00;
                              }
        /* other links in the sidebarmenu */
        .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                              background-color: #001787;
                              color: #ffffff;
        }
        /* other links in the sidebarmenu when hovered */
         .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                              background-color: #8d89ad;
         }
        /* toggle button when hovered  */
         .skin-blue .main-header .navbar .sidebar-toggle{
                              color:#c4ff00;
                              }
        /* toggle button when hovered  */
         .skin-blue .box.box-solid.box-primary>.box-header {
                              background-color: #001787;
                              color:#c4ff00;
         }

                              '))),

    autoWaiter(c("n_histogramm","n_distribution","n_exampletable","xmltodata")),
    tabItems(
      #Home----
           tabItem(tabName="home",
                fluidRow(
                  column(width=4,
                  box(width=12,
                    h1("MoonNMR - Local Edition"),
                    #img(src="https://user-images.githubusercontent.com/88379260/157672281-8f3902d3-998e-48cc-a445-25dc17a42fa5.png", height="20%",width="20%",align="left"),
                    h3(HTML("<b>M</b>etabol<b>O</b>mics <b>O</b>rga<b>N</b>iser NMR"),align="left"),
                    p("Tools for high-throughput NMR and untargeted Metabolomics",align="left"),
                    br(),
                    p("This is the local version of MoonNMR for the online version",a("click here",href="https://funkam.shinyapps.io/MoonShiny/")),
                    br(),
                    p("Developed by",a("Alexander Funk",href="https://www.uniklinikum-dresden.de/de/das-klinikum/kliniken-polikliniken-institute/klinische-chemie-und-laboratoriumsmedizin/forschung/copy_of_EMS")),
                    p("Institute for Clinical Chemistry and Laboratory Medicine"),
                    p("University Hospital Dresden, Fetscherstr. 74, 01307 Dresden")
                  ),
                  tabBox(width=12,title="Tools",
                         tabPanel(title="Preparation",
                                  h4(HTML("<b>IVDr Template Creator</b>")),
                                  p("Create ICON template for use with IVDr protocols and Bruker SampleJet, either by upload or by sample-by-sample creation."),
                                  p("Includes automatic filling of sample spaces in racks (1-96), automatic change of rack spaces (1-5) and interactive experiment selection."),
                                  p("Allows quick creation of large sample tables ready for submission."),
                                  h4(HTML("<b>ICON Template Creator</b>")),
                                  p("An additinal Tool for ICON is available that allows non-IVDr template creation with an upload of a list of experiments to be performed."),
                         ),
                         tabPanel(title="Processing",
                                  h4(HTML("<b>Import .xml into csv</b>")),
                                  p("Upload a zipped folder containing multiple .xml from B.I.QUANT-Methods and B.I.LISA"),
                                  p("Create an organised and cleaned-up datatable"),
                                  h4(HTML("<b>Clean-Up machine-CSV</b>")),
                                  p("Extract data from Bruker created files. Removes unwanted columns and uses HMBD names, when possible."),
                                  p("Due to the uneven text output by Bruker of Metabolic Profiles (;-csv) and Lipidomics Analysis (Tab-csv), it is possible to upload either one"),
                                  h4(HTML("<b>Combine tables</b>")),
                                  p("Combine two tables when a shared column (UniqueID) is present. The columns can be selected from list."),
                         ),
                         tabPanel(title="Analysis",
                                  h4(HTML("<b>Normalization</b>")),
                                  p("Normalize data in preparation for parametrric statistical analysis. Histogramm and Distrubtion plot interactively change upon option selection."),
                                  h4(HTML("<b>Automated Bar Plots</b>")),
                                  p("Upload a data table and interactively select columns to plot against a group column. Creation of large amount of plots possible, but are stored in a single PDF."),
                                  p("Suggested grouped plots (e.g. Lipid-Subclasses involving Triglycerides) are also available.")
                         )
                  ),
                  # box(width=12,
                  #   title="Changelog",status="primary",solidHeader=TRUE,
                  #   tableOutput("changelog")
                  # )
                  ),
                  column(width=8,
                        tabBox(width=12,title="Archive Plots",selected="Samples",
                          tabPanel(title="Samples",
                                   plotlyOutput(width="100%","archive_2020")),
                          tabPanel(title="Projects",
                                   plotlyOutput(width="100%","archive_project_2020"))
                        ),
                        box(width=12,title="Archive",status="primary",solidHeader = TRUE,
                            DT::dataTableOutput('archive')
                        )
                  )
                )

           ),

      #IVDr Uploader----
      tabItem(tabName="creator",
              h2("Create IVDr ICON Template from File"),
              fluidRow(
                box(width=3,title="Input",status="primary",solidHeader=TRUE,
                    radioGroupButtons("file1_type_Input","Choose File type",
                                      choices = list(".csv/txt" = 1, ".xlsx" = 2),
                                      selected = 1),
                    fileInput("file1", "File input", multiple=FALSE),
                    #textInput('file1sheet','Name of Sheet (Case-Sensitive)'),
                    radioGroupButtons("type","Sample Type",
                                      choices=c("Plasma"="plasma","Urine(neo)"="urine_neo","Urine"="urine","CSF"="CSF","MEOH"="MEOH"),
                                      selected="plasma"),
                    radioGroupButtons("size","Sample Size",
                                      choices=c("5 mm"="5mm","3 mm"="3mm"),
                                      selected="5mm"
                    ),
                    pickerInput("experiment", "Select experiments",
                                choices=c(),
                                options=list("actions-box"=TRUE, size=15,`selected-text-format` = "count > 3"),multiple=TRUE,selected=c(),
                    ),
                    radioGroupButtons("rack","Prefered Rack",
                                      choices=c("1"=1,"2"=2,"3"=3,"4"=4,"5"=5),
                                      selected="1"),
                    textInput("path", "Path",value="D:/data"),
                    selectizeInput("project", "Project", choices = startProjects, options=list(create=TRUE)),

                    tags$hr(),
                    actionButton("archive", "Add to Archive"),br(),br(),
                    downloadButton("downloadData", "Download Template")

                ),
                tabBox(width=8,selected="Template",title="Output",
                       tabPanel("Template", tableOutput("template")),
                       tabPanel("Raw Data", tableOutput("raw"))
                )
              )
      ),
      #IVDr Manual----
          tabItem(tabName="manual",
                  h2("Create manual (sample by sample) IVDr template for ICON"),
                  fluidRow(
                    box(width=3,title="Input",status="primary", solidHeader=TRUE,
                      radioGroupButtons("m_type","Sample Type",
                                   choices=c("Plasma"="m_plasma","Urine"="m_urine","Urine(neo)"="m_urine_neo","CSF"="m_csf","MEOH"="m_meoh"),
                                   selected="m_plasma"),
                      radioGroupButtons("m_size","Sample Size",
                                   choices=c("5 mm"="m_5mm","3 mm"="m_3mm"),
                                   selected="m_5mm"),
                      pickerInput("m_experiment", "Select experiments",
                                  choices=c(),
                                  options=list("actions-box"=TRUE, size=15,`selected-text-format` = "count > 3"),multiple=TRUE
                      ),
                      radioGroupButtons("m_rack","Prefered Rack",
                                   choices=c(1,2,3,4,5),
                                   selected=1),
                      numericInput("m_slot", "Position", value=1, min=1, max=96, step=1),
                      selectizeInput("m_project", "Project", choices = startProjects, options=list(create=TRUE)),
                      textInput("m_samplename","Sample Name"),
                      textInput("m_path", "Path",value="D:/IVDrData/data"),
                      tags$hr(),
                      actionButton("m_add", "Add Sample to Table"),
                      tags$hr(),
                      downloadButton("m_downloadData", "Download Template")
                    ),
                    box(width=7,title="Output",status="primary",solidHeader=TRUE, tableOutput("m_table2")
                    ),
                    )

          ),
      # #ICON Uploader----
      tabItem(tabName="iconup",
              h2("Create ICON Template from File"),
              fluidRow(
                box(width=3,title="Input",status="primary",solidHeader=TRUE,
                    radioGroupButtons("i_file1_type_Input","Choose File type",
                                      choices = list(".csv/txt" = 1, ".xlsx" = 2),
                                      selected = 1),
                    fileInput("i_file1", "File input", multiple=FALSE),
                    #textInput('file1sheet','Name of Sheet (Case-Sensitive)'),
                    textInput("i_type","Solvent"),
                    radioGroupButtons("i_size","Sample Size",
                                      choices=c("5 mm"="5mm","3 mm"="3mm"),
                                      selected="5mm"
                    ),
                    radioGroupButtons("i_file2_type_Input","Choose File type for Experiment List",
                                      choices = list(".csv/txt" = 1, ".xlsx" = 2),
                                      selected = 1),
                    fileInput("i_file2", "Upload Experiment List (single column file)", multiple=FALSE),
                    pickerInput("i_experiment", "Select experiments",
                                choices=c(),
                                options=list("actions-box"=TRUE, size=15,`selected-text-format` = "count > 3"),multiple=TRUE
                    ),
                    radioGroupButtons("i_rack","Prefered Rack",
                                      choices=c("1"=1,"2"=2,"3"=3,"4"=4,"5"=5),
                                      selected="1"),
                    textInput("i_path", "Path",value="D:/data"),
                    selectizeInput("i_project", "Project", choices = startProjects, options=list(create=TRUE)),

                    tags$hr(),
                    actionButton("i_archive", "Add to Archive"),br(),br(),
                    downloadButton("i_downloadData", "Download Template")

                ),
                tabBox(width=8,selected="Template",title="Output",
                       tabPanel("Template", tableOutput("i_template")),
                       tabPanel("Raw Data", tableOutput("i_raw"))



                )
              )
      ),
      # #ICON Creator----
      tabItem(tabName="iconman",
              h2("Create manual (sample by sample) IVDr template for ICON"),
              fluidRow(
                box(width=3,title="Input",status="primary", solidHeader=TRUE,
                    textInput("im_type","Solvent"),
                    radioGroupButtons("im_size","Sample Size",
                                      choices=c("5 mm"="5mm","3 mm"="3mm"),
                                      selected="5mm"),
                    radioGroupButtons("im_file2_type_Input","Choose File type for Experiment List",
                                      choices = list(".csv/txt" = 1, ".xlsx" = 2),
                                      selected = 1),
                    fileInput("im_file2", "Upload Experiment List (single column file)", multiple=FALSE),
                    pickerInput("im_experiment", "Select experiments",
                                choices=c(),
                                options=list("actions-box"=TRUE, size=15,`selected-text-format` = "count > 3"),multiple=TRUE
                    ),
                    radioGroupButtons("im_rack","Prefered Rack",
                                      choices=c(1,2,3,4,5),
                                      selected=1),
                    numericInput("im_slot", "Position", value=1, min=1, max=96, step=1),
                    selectizeInput("im_project", "Project", choices = startProjects, options=list(create=TRUE)),
                    textInput("im_samplename","Sample Name"),
                    textInput("im_path", "Path",value="D:/data"),
                    tags$hr(),
                    actionButton("im_add", "Add Sample to Table"),
                    tags$hr(),
                    downloadButton("im_downloadData", "Download Template")
                ),
                box(width=7,title="Output",status="primary",solidHeader=TRUE, tableOutput("im_table2")
                ),
              )

      ),
      #xml-er----
      tabItem(tabName="xml",
              h2("Extract parameters directly from zipped .xml files"),
              fluidRow(
                box(width=2,title="Input",status="primary",solidHeader=TRUE,
                    radioGroupButtons("xml_type","Choose panel",
                                      choices=c("Metabolites"="xml_metas","Lipids"="xml_lipids"),
                                      selected="xml_metas"),
                    fileInput("xml_file", "Choose zipped file", multiple=FALSE)
                ),
                box(width = 2,title="Options",status="primary",solidHeader=TRUE,

                    pickerInput("xml_keepers", "Select columns to keep",
                                choices=c("_conc","_errConc","_lod","rawConc","sigCorr"),
                                options=list("actions-box"=TRUE, size=15,`selected-text-format` = "count > 3"
                                ),
                                multiple=TRUE,selected=c("_conc","_errConc","_lod","rawConc","sigCorr")),
                    p("Only applicaple to Metabolite Panels"),
                ),

                box(width=2,title="Download",status="primary",solidHeader=TRUE,
                    downloadButton("xml_downloadData","Download table as .csv")
                ),
                box(width=4, background="black",color="yellow",
                    p("Upload a .zip file with multiple .xml reports of either Metabolites (Plasma/Urine) or a Lipid Report.",style="color:yellow"),
                    p("File names of the indivdiaul files are irrelevant.",style="color:yellow"),
                    p("A data table is created and can be downloaded.",style="color:yellow"),
                    p("Depending on the amount of xml files, this can take a few seconds.",style="color:yellow")

                )
              ),
              fluidRow(
                box(width=12, title="Table",status="primary",solidHeader=TRUE,
                    div(style = 'overflow-x: scroll',  DT::dataTableOutput('xmlfinal'))
                )
              )
      ),


      #Extractor----
          tabItem(tabName="extractor",
                  h2("Extract data from Bruker outputs"),
                  fluidRow(
                    box(width=3,title="Input",status="primary",solidHeader=TRUE,
                      radioGroupButtons("extractor_type","CSV Typ",
                                     choices=c(";"="e_colon",","="e_comma","Tab"="e_tab")),
                      p("Metabolite: Semicolon"),
                      p("Lipoproteine: Tab"),
                      radioGroupButtons("extractor_names", "Use standard names?",
                                   choices=c("Yes"="extractor_yes","No"="extractor_no")
                                   ),
                      fileInput("extractor_file", "Pick file", multiple=FALSE),
                      hr(),
                      downloadButton("downloadcsv", "Download CSV")
                    ),
                    box(width=9,title="Output",selected="Metabolite",status="primary",solidHeader=TRUE,
                      tabPanel("Metabolite",
                           div(style = 'overflow-x: scroll',  DT::dataTableOutput('final_csv'))
                              )
                    )
                  ),
                  fluidRow(
                  box(width=3, background="black",color="yellow",
                      p("Clean-Up the extracted text files from Bruker's IVDr Data Browser",style="color:yellow"),
                      p("Names can be standardized, sample names are cleaned up.",style="color:yellow")
                  )
                  )
          ),
      #deleter----
      tabItem(tabName = "deleter",
              fluidRow(
                box(width=3,title="Upload data:",status="primary",solidHeader=TRUE,
                    radioGroupButtons("del_file_type_Input","Choose File type",
                                      choices = list(".csv/txt" = 1, ".xlsx" = 2),
                                      selected = 1
                    ),
                    fileInput("del_file", "File input", multiple=FALSE),
                ),
                box(width=4,title="Options:",status="primary",solidHeader=TRUE,
                    p("Select type:"),
                    tabBox(width=12,selected="Parameter",
                           tabPanel("Parameter",
                                    pickerInput("del_column", "Choose parameters (columns) to delete",
                                                choices=c(),
                                                options=list("actions-box"=TRUE, size=15,`selected-text-format` = "count > 3"
                                                ),multiple=TRUE
                                    )
                           ),
                           tabPanel("Sample",
                                    pickerInput("del_irow", "Choose identifier column (e.g. SampleID)",
                                                choices=c(""),
                                                options=list("actions-box"=TRUE, size=15,`selected-text-format` = "count > 3"
                                                ),multiple=FALSE
                                    ),
                                    pickerInput("del_row", "Choose samples (rows) to delete",
                                                choices=c(),
                                                options=list("actions-box"=TRUE, size=15,`selected-text-format` = "count > 3"
                                                ),multiple=TRUE
                                    )
                                    
                           ),
                           tabPanel("Group",
                                    pickerInput("del_igroup", "Choose identifier column (e.g. Group)",
                                                choices=c(),
                                                options=list("actions-box"=TRUE, size=15,`selected-text-format` = "count > 3"
                                                ),multiple=FALSE
                                    ),
                                    pickerInput("del_group", "Choose groups to delete",
                                                choices=c(),
                                                options=list("actions-box"=TRUE, size=15,`selected-text-format` = "count > 3"
                                                ),multiple=TRUE
                                    )
                           )
                    )
                ),
                box(width=2,title="Download:",status="primary",solidHeader=TRUE,
                    downloadButton("del_downloader","Download table"),
                ),
                box(width=5, background="black",
                    p("Upload a data table and edit the parameters, samples and groups interactively.",style="color:yellow"),
                    p("Multiple edits possible at the same time.",style="color:yellow")
                    
                )
              ),
              
              fluidRow(
                box(width=12, title="Output",status="primary",solidHeader=TRUE,
                    div(style = 'overflow-x: scroll',  DT::dataTableOutput('del_table'))
                )
              )
              
      ),
      
      
      
      
      #Combiner-----
      tabItem(tabName="combiner",
              fluidRow(
                box(width=2,title="Upload first Table:",status="primary",solidHeader=TRUE,
                    radioGroupButtons("c_fileleft_type_Input","Choose file type",
                                      choices = list(".csv/txt" = 1, ".xlsx" = 2),
                                      selected = 1
                    ),
                    fileInput("c_fileleft", "File input", multiple=FALSE),
                    pickerInput("c_columnleft", "Select column of first table",
                                choices=c(),
                                options=list("actions-box"=TRUE, size=15,`selected-text-format` = "count > 3"
                                ),
                                multiple=FALSE),
                ),

                box(width=2,title="Upload second Table:",status="primary",solidHeader=TRUE,
                    radioGroupButtons("c_fileright_type_Input","Choose file type",
                                      choices = list(".csv/txt" = 1, ".xlsx" = 2),
                                      selected = 1
                    ),
                    fileInput("c_fileright", "File input", multiple=FALSE),
                    pickerInput("c_columnright", "Select column of second table",
                                choices=c(),
                                options=list("actions-box"=TRUE, size=15,`selected-text-format` = "count > 3"
                                ),
                                multiple=FALSE)
                ),
                box(width=2,title="Download:",status="primary",solidHeader=TRUE,
                    useWaiter(),
                    downloadButton("c_downloader","Download joined table")
                ),
                box(width=3, background="black",
                    p("Combine two tables with an identical column, ideally a SampleID",style="color:yellow"),
                    p("The columns do not need to have the same name",style="color:yellow"),
                    p("The column in the output table will have the name UniqueID",style="color:yellow")
                )
              ),
              fluidRow(
                tabBox(width=12,
                       tabPanel("Table 1",
                                div(style = 'overflow-x: scroll',  DT::dataTableOutput('c_tableleft'))
                       ),
                       tabPanel("Table 2",
                                div(style = 'overflow-x: scroll',  DT::dataTableOutput('c_tableright'))
                       ),
                       tabPanel("Table combined",
                                div(style = 'overflow-x: scroll',  DT::dataTableOutput('c_tablecombined'))
                       )
                )
              )
      ),
      #Normalizor----
         tabItem(tabName="normalizer",
                 h2("Normalize data"),

            fluidRow(
                   box(width=2,title="Input",status="primary",solidHeader=TRUE,
                      radioGroupButtons("n_file_type_Input","Choose File type",
                                        choices = list(".csv/txt" = 1, ".xlsx" = 2),
                                        selected = 1
                      ),
                      fileInput("n_file", "File input", multiple=FALSE),
                      radioGroupButtons("n_log","Log10",
                                        choices=c("Yes"="n_log_yes","No"="n_log_no"),
                                        selected="n_log_no"
                      ),
                      radioGroupButtons("n_center","MeanCenter",
                                        choices=c("Yes"="n_center_yes","No"="n_center_no"),
                                        selected="n_center_no"
                      ),
                      downloadButton("n_download", "Download CSV")
                   ),
                  box(width=5,title="Histogramm",status="primary",solidHeader=TRUE,
                     plotlyOutput(width="100%","n_histogramm")
                  ),
                  box(width=5,title="Distribution",status="primary",solidHeader=TRUE,
                     div(style="height:400px;overflow-y:scroll",plotlyOutput(width="100%", height="2000px","n_distribution"))
                  )
            ),
            fluidRow(
                   box(width=12,title="Output",status="primary",solidHeader=TRUE,
                      div(style = 'overflow-x: scroll',  DT::dataTableOutput('n_exampletable'))
                   )
            )
       ),
      #Plotting----
      tabItem(tabName="plots",
              h2("Create plots of column vs a group/time"),
              fluidRow(
                column(width=3,
                       box(width=12,title="Upload data:",status="primary",solidHeader=TRUE,
                           radioGroupButtons("p_file_type_Input","Choose File type",
                                             choices = list(".csv/txt" = 1, ".xlsx" = 2),
                                             selected = 1
                           ),
                           fileInput("p_file", "File input", multiple=FALSE),
                           hr(),
                           pickerInput("p_group", "Select column for grouping",
                                       choices=c(),
                                       options=list("actions-box"=TRUE, size=15,`selected-text-format` = "count > 3"
                                       ),
                                       multiple=FALSE),
                           pickerInput("p_column", "Choose columns to plot",
                                       choices=c(),
                                       options=list("actions-box"=TRUE, size=15,`selected-text-format` = "count > 3"
                                       ),
                                       multiple=TRUE),
                           hr(),
                           shinyjqui::orderInput("p_groupsorter","Drag to re-order groups/timepoints",
                                      items=c())
                       )
                ),
                column(width=4,
                       fluidRow(
                         box(width=12,title="Chose plot Options:",status="primary",solidHeader=TRUE,
                             radioGroupButtons("p_stat","Statistics",
                                               choices=c("None"="p_none","t-test"="p_ttest","Wilcoxon"="p_wilc","Anova"="p_anova","Kruskal-Wallis"="p_krusk"),
                                               selected="p_none",individual=TRUE)
                         )
                       ),
                       fluidRow(
                         box(width=12, background="black",
                             p("Individual plots are shown on the right, only one at a time.",style="color:yellow"),
                             p("Statistical significance is highlighted in the plots by *.",style="color:yellow"),
                             p("Interactive single plot on the right does not show significance.",style="color:yellow")
                             
                         )
                       ),
                       fluidRow(
                         box(width=12,title="Download:",status="primary",solidHeader=TRUE,
                             useWaiter(),
                             h5("Download a single or multiple individual plots as .PDF or a small selection of different parameters in a grouped plot (as .tiff)"),
                             hr(),
                             downloadButton("p_downloader","PDF"),
                             downloadButton("p_downloader_groups","Grouped")
                         )
                       )
                ),
                column(width=5,
                       box(width=12, status="primary",solidHeader=TRUE,title="Single plot display",
                           plotlyOutput("p_single"),
                           textOutput("p_text"),
                           textOutput("p_text2")
                       )
                )
              ),
              fluidRow(
                box(width=12, status="primary",solidHeader=TRUE,title="Data Table Overview",
                    div(style = 'overflow-x: scroll',DT::dataTableOutput("p_exampletable"))
                )
              )
      ),
      
    #PCA----
      tabItem(tabName="pca",
              h2("Principal Component Analysis"),
              fluidRow(
                box(width=3,title="Upload data:",status="primary",solidHeader=TRUE,
                    radioGroupButtons("pca_file_type_Input","Choose File type",
                                      choices = list(".csv/txt" = 1, ".xlsx" = 2),
                                      selected = 1
                    ),
                    fileInput("pca_file", "File input", multiple=FALSE),
                    pickerInput("pca_id", "Select Sample ID column",
                                choices=c(),
                                options=list("actions-box"=TRUE, size=15,`selected-text-format` = "count > 3"
                                ),
                                multiple=FALSE,selected=c()),
                    pickerInput("pca_group", "Select Group column",
                                choices=c(),
                                options=list("actions-box"=TRUE, size=15,`selected-text-format` = "count > 3"
                                ),
                                multiple=FALSE,selected=NULL),
                    radioGroupButtons("pca_labels","Include labels in plot",
                                      choices=c("Yes"="pca_labels_yes","No"="pca_labels_no"),
                                      selected="pca_labels_no"),
                    hr(),
                    downloadButton("pca_download","Download PCA table")
                ),
                tabBox(width=9,selected="Data Table",
                       tabPanel("Data Table",
                                div(style = 'overflow-x: scroll',DT::dataTableOutput("pca_table"))
                       ),
                       tabPanel("PCA",
                                div(style = 'overflow-x: scroll',DT::dataTableOutput("pca_pca"))
                       )
                )
              ),
              fluidRow(
                box(width=6, status="primary",solidHeader=TRUE,title="PCA Scores",
                    plotlyOutput(width="90%","pca_scores")
                ),
                box(width=6, status="primary",solidHeader=TRUE,title="PCA Biplot",
                    plotlyOutput(width="90%","pca_biplot")
                )
              )
      )
  )
  )
)






# Server----
app_server <- function(input,output,session) {

#Server Home----
  #create reactive archive, ->in excel tabs
  archive<-reactive({read.csv("archive.csv")})
  #Output
  output$archive<-DT::renderDataTable(archive())
  #output archive plots
  output$archive_2020 <- renderPlotly({
    archive4plot<-archive()
    archive4plot$Date<-format(as.Date(archive4plot$Date, format="%m/%d/%Y"))
    ggplot(archive4plot, aes(x=Type, fill=Type))+geom_bar()+theme_few()+theme(legend.position="none",axis.title.y=element_blank())+labs(x="Sample",y="Samples")+coord_flip()})
  output$archive_project_2020 <- renderPlotly({
    archive4plot<-archive()
    archive4plot$Date<-format(as.Date(archive4plot$Date, format="%m/%d/%Y"))
    ggplot(archive4plot, aes(x=Project, fill=Project))+geom_bar()+theme_few()+theme(legend.position="none",axis.title.y=element_blank())+labs(y="Samples")+coord_flip()})

#Server Excel Importer----
  solvent<-reactive({switch(input$type, "urine_neo"="Urine_Neo","urine"="Urine", "plasma"="Plasma","CSF"="CSF","MEOH"="MEOH")})
  rack<-reactive({switch(input$rack, "1"=1,"2"=2,"3"=3,"4"=4,"5"=5)})
  size<-reactive({switch(input$size, "5mm"="5mm","3mm"="3mm")})
  path<-reactive(input$path)
  project<-reactive({input$project})


  raw_data<-reactive({
    inFile <- input$file1
    if (is.null(inFile)) {
      return(NULL) }
    if (input$file1_type_Input == "1") {
      read.csv(inFile$datapath,
               header = TRUE,
               stringsAsFactors = FALSE)
    } else {
      read.xlsx(inFile$datapath)
    }
  })

  experiments<-reactive({input$experiment})

  elementpicker<-reactive({
    solvent<-solvent()
    size<-size()
    elementpicked<-experimentpicker(solvent,size)
  })

  #output$elementtext<-renderText({elementpicker()})

  #update picks
  experiment<-observe({
    req(input$file1)
    updatePickerInput(
      session,
      "experiment",
      choices = elementpicker()
    )
  })

  #Create reactive template
  create_template<-reactive({
    data_output<-excel_creator(raw_data(),solvent(),size(),rack(),experiments(),path())
  })
  #create archive addition
  archive_add<-reactive({
    Name<-raw_data()
    archive_current<-archive()
    Project<-project()
    Date<-Sys.Date()
    Type<-solvent()
    Size<-size()
    Date<-as.character(Date) #likely needs adjusting
    archive_updated<-cbind(Date,Name,Project,Type,Size)
    new_archive<-rbind(archive_current,archive_updated)
  })
  observeEvent(input$archive, {
    write.csv(archive(),"archive_backup.csv",row.names=FALSE)
    write.csv(archive_add(),"archive.csv",row.names=FALSE)
    showModal(modalDialog(
      title = "Important message",
      "The samples have been added to the archive."
    ))
  })


      #Output
          output$template<-renderTable({
            req(input$file1)
            create_template()
          })
          output$raw<-renderTable({
            req(input$file1)
            raw_data()
          })
          output$changelog<-renderTable({
            changelog
          })

      #Buttons
          output$downloadData <- downloadHandler(
            filename = function() {
              paste("Template_", ".xlsx", sep = "")
            },
            content = function(file) {
              write.xlsx(create_template(), file, row.names = FALSE)
            }
          )


#Server Excel Manual----
        #reactives
          m_solvent<-reactive({switch(input$m_type, "m_urine"="Urine", "m_plasma"="Plasma","m_media"="Media","m_urine_neo"="Urine_Neo","m_csf"="CSF","m_meoh"="MEOH")})
          m_rack<-reactive({switch(input$m_rack, "1"=1,"2"=2,"3"=3,"4"=4,"5"=5)})
          m_size<-reactive({switch(input$m_size, "m_5mm"="5mm","m_3mm"="3mm")})
          m_samplename<-reactive({input$m_samplename})
          m_slot<-reactive({input$m_slot})
          m_path<-reactive({input$m_path})
          m_experiments<-reactive({input$m_experiment})
          m_project<-reactive({input$m_project})
          m_archive<-reactive({read.csv("archive.csv")})


          m_elementpicker<-reactive({
            solvent<-m_solvent()
            size<-m_size()
            elementpicked<-experimentpicker(solvent,size)
          })

          #output$elementtext<-renderText({elementpicker()})

          #update picks
          m_experiment<-observe({
            req(input$m_type)
            updatePickerInput(
              session,
              "m_experiment",
              choices = m_elementpicker()
            )
          })
          manual_archive<-reactive({
            Date<-as.character(Sys.Date())
            Name<-m_samplename()
            Project<-m_project()
            Type<-m_solvent()
            Size<-m_size()
            m_archiver<-cbind(Date,Name,Project,Type,Size)
          })
          m_updater_archive<-eventReactive(input$m_add,{
            m_looper_archive<<-rbind(m_looper_archive,manual_archive())
          })
        #read in manual empty template
          m_looper<-manual_template<-data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("DISK","NAME","SOLVENT","EXPERIMENT","HOLDER", "TITLE"))))
          m_looper_archive<-data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("Date","Name","Project","Type","Size"))))
        #create_manual_method
          manual_samples<-reactive({
            m_new_sample<-m_samplename()
            m_new_sample<-as.data.frame(m_new_sample)
            names(m_new_sample)[names(m_new_sample) =="m_new_sample"]<-"Name"
            m_table<-excel_manual(m_new_sample,m_solvent(),m_size(),m_rack(),m_slot(),m_experiments(),m_path())
          })
        #append dataframes
          m_updater <- eventReactive(input$m_add,{
            m_slot_counter <- m_slot()+1
            if(m_slot_counter>96){
              m_slot_counter<-1
              updateNumericInput(session, "m_slot", value = m_slot_counter)
              m_rack_counter<-m_rack()+1
              if(m_rack_counter>5){
                m_rack_counter<-1
                updateRadioButtons(session, "m_rack", choices=c(1,2,3,4,5), selected=m_rack_counter)
              } else {
                updateRadioButtons(session, "m_rack", choices=c(1,2,3,4,5), selected=m_rack_counter)
              }
            } else {
              updateNumericInput(session, "m_slot", value = m_slot_counter)
            }
            write.csv(m_archive(),"archive_backup.csv",row.names=FALSE)
            write.csv(m_archive_add(),"archive.csv",row.names=FALSE)
            m_looper<<-rbind(m_looper,manual_samples())
          })
          m_archive_add<-reactive({
            archive_current<-archive()
            new_archive<-rbind(archive_current,m_updater_archive())
          })
        #Output
          output$m_table2<-renderTable({m_updater()})
        #Buttons
             output$m_downloadData <- downloadHandler(
                filename = function() {
                paste("Template_manual", ".xlsx", sep = "")
              },
              content = function(file) {
                write.xlsx(m_updater(), file, row.names = FALSE)
              }
            )
             observeEvent(input$m_archive, {
               write.csv(m_archive(),"archive_backup.csv",row.names=FALSE)
               write.csv(m_archive_add(),"archive.csv",row.names=FALSE)
               showModal(modalDialog(
                 title = "Important message",
                 "The samples have been added to the archive."
               ))
             })

#Server ICON Uploader----
             i_solvent<-reactive({input$i_type})
             i_rack<-reactive({switch(input$i_rack, "1"=1,"2"=2,"3"=3,"4"=4,"5"=5)})
             i_size<-reactive({switch(input$i_size, "5mm"="5mm","3mm"="3mm")})
             i_path<-reactive({input$i_path})
             i_experiments<-reactive({input$i_experiment})
             i_project<-reactive({input$project})
             i_archive<-reactive({read.csv("archive.csv")})



             i_raw_data<-reactive({
               inFile <- input$i_file1
               if (is.null(inFile)) {
                 return(NULL) }
               if (input$i_file1_type_Input == "1") {
                 read.csv(inFile$datapath,
                          header = TRUE,
                          stringsAsFactors = FALSE)
               } else {
                 read.xlsx(inFile$datapath)
               }
             })

             i_experimentlist<-reactive({
               inFile <- input$i_file2
               if (is.null(inFile)) {
                 return(NULL) }
               if (input$i_file2_type_Input == "1") {
                 read.csv(inFile$datapath,
                          header = TRUE,
                          stringsAsFactors = FALSE)
               } else {
                 read.xlsx(inFile$datapath)
               }
             })


             #update picks
             i_experiment<-observe({
               req(input$i_file2)
               updatePickerInput(
                 session,
                 "i_experiment",
                 choices = i_experimentlist()
               )
             })

             #Create reactive template
             i_create_template<-reactive({
               data_output<-excel_creator(i_raw_data(),i_solvent(),i_size(),i_rack(),i_experiments(),i_path())
             })
             #Output
             output$i_template<-renderTable({
               req(input$i_file1,input$i_file2)
               i_create_template()
             })
             output$i_raw<-renderTable({
               req(input$i_file1)
               i_raw_data()
             })

             #create archive addition
             i_archive_add<-reactive({
               Name<-i_raw_data()
               archive_current<-i_archive()
               Project<-i_project()
               Date<-Sys.Date()
               Type<-i_solvent()
               Size<-i_size()
               Date<-as.character(Date) #likely needs adjusting
               archive_updated<-cbind(Date,Name,Project,Type,Size)
               new_archive<-rbind(archive_current,archive_updated)
             })
             observeEvent(input$i_archive, {
               write.csv(i_archive(),"archive_backup.csv",row.names=FALSE)
               write.csv(i_archive_add(),"archive.csv",row.names=FALSE)
               showModal(modalDialog(
                 title = "Important message",
                 "The samples have been added to the archive."
               ))
             })


             #Buttons
             output$i_downloadData <- downloadHandler(
               filename = function() {
                 paste("Template_", ".xlsx", sep = "")
               },
               content = function(file) {
                 write.xlsx(i_create_template(), file, row.names = FALSE)
               }
             )


#Server ICON Manual----
             #reactives
             im_solvent<-reactive(input$im_type)
             im_rack<-reactive({switch(input$im_rack, "1"=1,"2"=2,"3"=3,"4"=4,"5"=5)})
             im_size<-reactive({switch(input$im_size, "5mm"="5mm","3mm"="3mm")})
             im_samplename<-reactive({input$im_samplename})
             im_slot<-reactive({input$im_slot})
             im_path<-reactive({input$im_path})
             im_experiments<-reactive({input$im_experiment})
             im_project<-reactive({input$im_project})
             im_archive<-reactive({read.csv("archive.csv")})

             im_experimentlist<-reactive({
               inFile <- input$im_file2
               if (is.null(inFile)) {
                 return(NULL) }
               if (input$im_file2_type_Input == "1") {
                 read.csv(inFile$datapath,
                          header = TRUE,
                          stringsAsFactors = FALSE)
               } else {
                 read.xlsx(inFile$datapath)
               }
             })

             #update picks
             im_experiment<-observe({
               req(input$im_file2)
               updatePickerInput(
                 session,
                 "im_experiment",
                 choices = im_experimentlist()
               )
             })

             im_manual_archive<-reactive({
               Date<-as.character(Sys.Date())
               Name<-im_samplename()
               Project<-im_project()
               Type<-im_solvent()
               Size<-im_size()
               im_archiver<-cbind(Date,Name,Project,Type,Size)
             })

             im_updater_archive<-eventReactive(input$im_add,{
               im_looper_archive<<-rbind(im_looper_archive,im_manual_archive())
             })


             #read in manual empty template
             im_looper<-manual_template<-data.frame(matrix(ncol=6,nrow=0, dimnames=list(NULL, c("DISK","NAME","SOLVENT","EXPERIMENT","HOLDER", "TITLE"))))
             im_looper_archive<-data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("Date","Name","Project","Type","Size"))))

             #create_manual_method
             im_manual_samples<-reactive({
               im_new_sample<-im_samplename()
               im_new_sample<-as.data.frame(im_new_sample)
               names(im_new_sample)[names(im_new_sample) =="im_new_sample"]<-"Name"
               m_table<-excel_manual(im_new_sample,im_solvent(),im_size(),im_rack(),im_slot(),im_experiments(),im_path())
             })
             #append dataframes
             im_updater <- eventReactive(input$im_add,{
               im_slot_counter <- im_slot()+1
               if(im_slot_counter>96){
                 im_slot_counter<-1
                 updateNumericInput(session, "im_slot", value = im_slot_counter)
                 im_rack_counter<-im_rack()+1
                 if(im_rack_counter>5){
                   im_rack_counter<-1
                   updateRadioButtons(session, "im_rack", choices=c(1,2,3,4,5), selected=im_rack_counter)
                 } else {
                   updateRadioButtons(session, "im_rack", choices=c(1,2,3,4,5), selected=im_rack_counter)
                 }
               } else {
                 updateNumericInput(session, "im_slot", value = im_slot_counter)
               }
               write.csv(im_archive(),"archive_backup.csv",row.names=FALSE)
               write.csv(im_archive_add(),"archive.csv",row.names=FALSE)
               im_looper<<-rbind(im_looper,im_manual_samples())
             })

             im_archive_add<-reactive({
               im_archive_current<-im_archive()
               new_archive<-rbind(im_archive_current,im_updater_archive())
             })


             #Output
             output$im_table2<-renderTable({im_updater()})
             #Buttons
             output$im_downloadData <- downloadHandler(
               filename = function() {
                 paste("Template_manual", ".xlsx", sep = "")
               },
               content = function(file) {
                 write.xlsx(im_updater(), file, row.names = FALSE)
               }
             )

#Server XML----
             raw_xmllist<-reactive({
               inFile <- input$xml_file
               if (is.null(inFile)) {
                 return(NULL) }
               lapply(unzip(inFile$datapath),read_xml)
             })
             xml_keepers<-reactive(input$xml_keepers)
             xml_type<-reactive({switch(input$xml_type, "xml_metas"="Metabolites","xml_lipids"="Lipids")})

             xmlextractor<-reactive({
               xmlextractor<-xmler(raw_xmllist(),xml_type())
             })
             xmlcleaner<-reactive({
               if(xml_type()=="Metabolites"){
                 keepers<-xml_keepers()
                 df1<-xmlextractor()
                 Sample<-df1$Sample
                 combinekeepers<-as.data.frame(lapply(keepers, FUN=function(x) {
                   df1<-select(df1,ends_with(x))
                 }))
                 combinekeepers<-combinekeepers[,order(names(combinekeepers))]
                 combined<-cbind(Sample,combinekeepers)
               }
               else {df1<-xmlextractor()}
             })

             #Output
             output$xmlfinal<-DT::renderDataTable({
               req(input$xml_file)
               xmlcleaner()
             })

             output$xml_downloadData <- downloadHandler(
               filename = function() {
                 paste("Measurements_", ".csv", sep = "")
               },
               content = function(file) {
                 write.csv(xmlcleaner(), file, row.names = FALSE)
               }
             )

#Server Data Exctractor Tab----
        e_solvent<-reactive({switch(input$extractor_type, "e_comma"="Comma","e_colon"="Semicolon","e_tab"="Tab")})
        e_names<-reactive({switch(input$extractor_names, "extractor_yes"="Yes","extractor_no"="No")})
        dirty_csv<-reactive({

           if(e_solvent()=="Semicolon"){
            read.csv2(input$extractor_file$datapath)
           } else if(e_solvent()=="Comma"){
             read.csv(input$extractor_file$datapath)
           } else if(e_solvent()=="Tab"){
             read.csv(input$extractor_file$datapath,sep="\t")
           }

        })
        #create clean_csv
          cleaned_csv<-reactive({
            clean<-clean_csv(e_names(),dirty_csv())
            clean<-clean
          })
        #Output
          output$final_csv<-DT::renderDataTable({
            req(input$extractor_file)
            cleaned_csv()
          })
        #Buttons
          output$downloadcsv <- downloadHandler(
            filename = function() {
              paste("Cleaned_", ".csv", sep = "")
            },
            content = function(file) {
              write.csv(cleaned_csv(), file, row.names = FALSE)
            }
          )
#Server Deleter----
          del_file <- reactive({
            del_inFile <- input$del_file
            if (is.null(del_inFile))
              return(NULL)
            if (input$del_file_type_Input == "1") {
              read.csv(del_inFile$datapath,
                                   header = TRUE,
                                   stringsAsFactors = FALSE)
            } else {
              read.xlsx(del_inFile$datapath)
            }
          })
          
          del_columns<-reactive({input$del_column})
          del_rows<-reactive({input$del_row})
          del_irows<-reactive({input$del_irow})
          del_groups<-reactive({input$del_group})
          del_igroups<-reactive({input$del_igroup})
          samplerows<-reactive({
            req(del_irows())
            dat<-subset(del_file(),select=names(del_file()) %in% del_irows())
          })
          grouprows<-reactive({
            req(del_igroups())
            dat<-subset(del_file(),select=names(del_file()) %in% del_igroups())
            dat<-unique(dat)
          })
          #update picks
          
          del_column<-observe({
            updatePickerInput(
              session,
              "del_column",
              choices = names(del_file())
              
            )
          })
          del_irow<-observe({
            updatePickerInput(
              session,
              "del_irow",
              choices = names(del_file())
              
            )
          })
          del_row<-observe({
            req(input$del_file)
            req(del_irows())
            
            updatePickerInput(
              session,
              "del_row",
              choices = samplerows()
              
            )
          })
          del_igroup<-observe({
            updatePickerInput(
              session,
              "del_igroup",
              choices = names(del_file())
              
            )
          })
          del_group<-observe({
            req(input$del_file)
            req(del_igroups())
            updatePickerInput(
              session,
              "del_group",
              choices = grouprows()
            )
          })
          # del_table<-reactive({
          #   req(input$del_file)
          #
          #   if(length(del_columns()!=0)){
          #     "%ni%"<-Negate('%in%')
          #     del_final<-subset(del_file(),select=names(del_file()) %ni% del_columns())
          #   }
          
          
          # })
          del_table<-reactive({
            req(input$del_file)
            del_final<-del_file()
            if(sjmisc::is_empty(del_columns())==FALSE){
              "%ni%"<-Negate('%in%')
              del_final<-subset(del_final,select=names(del_final) %ni% del_columns())
            }
            if(sjmisc::is_empty(del_rows())==FALSE){
              irows<-as.name(del_irows())
              rows<-paste(del_rows(),collapse="|")
              del_final<-del_final %>% filter(!str_detect(!!as.symbol(irows),rows))
            }
            if(sjmisc::is_empty(del_groups())==FALSE){
              igroup<-as.name(del_igroups())
              group<-del_groups()
              del_final<-del_final %>% filter(!str_detect(!!as.symbol(igroup),group))
            }
            del_final<-del_final
          })
          
          #output
          output$del_table<-DT::renderDataTable(del_table())
          
          
          output$del_downloader <- downloadHandler(
            filename = function() {
              paste("Data_", ".csv", sep = "")
            },
            content = function(file) {
              write.csv(del_table(), file, row.names = FALSE)
            }
          )
          
#Server Combiner----
          c_fileleft <- reactive({
            inFilel <- input$c_fileleft
            if (is.null(inFilel))
              return(NULL)
            if (input$c_fileleft_type_Input == "1") {
              read.csv(inFilel$datapath,
                                   header = TRUE,
                                   stringsAsFactors = FALSE)
            } else {
              read.xlsx(inFile$datapath)
            }
          })
          c_fileright <- reactive({
            inFiler <- input$c_fileright
            if (is.null(inFiler))
              return(NULL)
            if (input$c_fileright_type_Input == "1") {
              read.csv(inFiler$datapath,
                                   header = TRUE,
                                   stringsAsFactors = FALSE)
            } else {
              read.xlsx(inFile$datapath)
            }
          })

          c_columnleft<-reactive({input$c_columnleft})
          c_columnright<-reactive({input$c_columnright})



          #update picks
          c_columnlefts<-observe({
            updatePickerInput(
              session,
              "c_columnleft",
              choices = names(c_fileleft())
            )
          })
          c_columnrights<-observe({
            updatePickerInput(
              session,
              "c_columnright",
              choices = names(c_fileright())
            )
          })

          c_tablejoined<-reactive({
            req(input$c_fileleft)
            req(input$c_fileright)

            dataleft<-c_fileleft()
            dataright<-c_fileright()
            groupleft<-c_columnleft()
            groupright<-c_columnright()
            group<-c("UniqueID")
            names(dataleft)[names(dataleft) == groupleft] <- group
            names(dataright)[names(dataright) == groupright] <- group
            data_merged<-merge(x=dataleft,y=dataright,by=group,all.x=TRUE)
          })

          #output
          output$c_tableleft<-DT::renderDataTable(c_fileleft())
          output$c_tableright<-DT::renderDataTable(c_fileright())
          output$c_tablecombined<-DT::renderDataTable(c_tablejoined())

          output$c_downloader <- downloadHandler(
            filename = function() {
              paste("Joined_csv", ".csv", sep = "")
            },
            content = function(file) {
              write.csv(c_tablejoined(), file, row.names = FALSE)
            }
          )


#Server Normalizer----
        n_log<-reactive({switch(input$n_log,"n_log_yes"="Yes","n_log_no"="No")})
        n_center<-reactive({switch(input$n_center,"n_center_yes"="Yes","n_center_no"="No")})

        n_file <- reactive({
          n_inFile <- input$n_file
          if (is.null(n_inFile))
            return(NULL)
          if (input$n_file_type_Input == "1") {
            read.csv(n_inFile$datapath,
                                 header = TRUE,
                                 stringsAsFactors = FALSE)
          } else {
            read.xlsx(n_inFile$datapath)
          }
        })
        w <- Waiter$new()

        n_normtable<-reactive({
          n_data_norm<-normalizor(n_file(),n_log(),n_center())

        })
        #Buttons
        output$n_download <- downloadHandler(
          filename = function() {
            paste("Norm_csv", ".csv", sep = "")
          },
          content = function(file) {
            write.csv(n_normtable(), file, row.names = FALSE)
          }
        )

        n_plotters<-reactive({
          req(input$n_file)
          df1<<-as.data.frame.array(melt(n_normtable()))
        })

        output$n_exampletable<-DT::renderDataTable({
          req(input$n_file)
          n_normtable()
        })
        output$n_histogramm<-renderPlotly({
          req(input$n_file)
          ggplot(n_plotters(),aes(x=value,fill=cut(value,100)))+geom_histogram()+theme_classic()+xlab("")+ylab("")+theme(legend.position = "none")
        })
        output$n_distribution<-renderPlotly({
          req(input$n_file)
          ggplot(n_plotters(),aes(x = value, y=variable, color=variable))+geom_line(size=3)+theme_classic()+theme(legend.position = "none",axis.text.y = element_text(size=10))+xlab("")+ylab("")

        })

#Server Plotting Tab----
        #input
        p_file <- reactive({
          inFile <- input$p_file
          if (is.null(inFile))
            return(NULL)
          if (input$p_file_type_Input == "1") {
            read.csv(inFile$datapath,
                                 header = TRUE,
                                 stringsAsFactors = FALSE)
          } else {
            read.xlsx(inFile$datapath)
          }
        })
        p_groupsorters<-reactive({input$p_groupsorter})
        p_columns<-reactive({input$p_column})
        p_groups<-reactive({input$p_group})
        p_grouped<-reactive({switch(input$p_grouped, "p_grouped_yes"="Yes","p_grouped_no"="No")})
        
        p_stat<-reactive({switch(input$p_stat,"p_none"="none","p_ttest"="ttest","p_wilc"="wilcox","p_anova"="anova","p_krusk"="kruskal")})
        
        p_uniquegroups<-reactive({
          req(p_file(),p_groups())
          pdat<-subset(p_file(),select=names(p_file()) %in% p_groups())
          pdat<-unique(pdat)
          pdat<-pdat[[p_groups()]]
        })
        #update picks
        p_column<-observe({
          updatePickerInput(
            session,
            "p_column",
            choices = names(p_file())
          )
        })
        p_group<-observe({
          updatePickerInput(
            session,
            "p_group",
            choices = names(p_file())
          )
        })
        p_groupsorter<-observe({
          shinyjqui::updateOrderInput(
            session,
            "p_groupsorter",
            items = p_uniquegroups()
          )
        })
        p_filesorted<-reactive({
          req(p_groups(),p_file(),p_groupsorters())
          p_groups<-as.symbol(p_groups())
          p_groupsorters<-p_groupsorters()
          p_file2<-p_file()
          p_file2<-mutate(p_file2, !!p_groups := !!p_groups %>% factor(levels = p_groupsorters))
        })
        p_singleplotter<-reactive({
          req(p_filesorted())
          plotter_single(p_filesorted(),p_groups(),p_columns(),p_stat())
        })

        
        #outputs
        output$p_exampletable<-DT::renderDataTable({
          req(input$p_file)
          p_file()
        })

        output$p_single <- renderPlotly({
          req(input$p_file,input$p_column)
          p_singleplotter()
        })
        output$p_downloader_groups <- downloadHandler(
          filename = "Grouped.tiff",
          content = function(file) {
            p_plot<-plotter_grouped(p_filesorted(),p_groups(),p_columns(),p_stat())
            tiff(file)
            print(p_plot)
            dev.off()
          })
        output$p_downloader <- downloadHandler(
          filename = "Plots_Summary.pdf",
          
          content = function(file) {
            waiter_show(html=waiting_screen,color=transparent(0.9))
            
            p_data<-p_filesorted()
            p_plot<-plotter_v2(p_data,p_groups(),p_columns(),p_stat())
            pdf(file)
            print(p_plot)
            dev.off()
            waiter_hide()
            sendSweetAlert(
              session = session,
              title = "Success!",
              text = "Plots generated!",
              type = "success"
            )
          })

#server PCA----
          pca_file <- reactive({
            inFile <- input$pca_file
            if (is.null(inFile))
              return(NULL)
            if (input$p_file_type_Input == "1") {
              read.csv(inFile$datapath,
                                   header = TRUE,
                                   stringsAsFactors = FALSE)
            } else {
              read.xlsx(inFile$datapath)
            }
          })
          pca_ids<-reactive({input$pca_id})
          pca_labels<-reactive({switch(input$pca_labels, "pca_labels_yes"="Yes","pca_labels_no"="No")})
          
          pca_groups<-reactive({input$pca_group})
          #update picks
          pca_group<-observe({
            updatePickerInput(
              session,
              "pca_group",
              choices = names(pca_file())
            )
          })
          pca_id<-observe({
            updatePickerInput(
              session,
              "pca_id",
              choices = names(pca_file())
            )
          })
          #pccomp
          pca_calc<-reactive({
            req(pca_ids())
            pca_data<-pca_file()
            #put sampleid as rownames for labels
            idloc<-grep(pca_ids(),colnames(pca_data))
            row.names(pca_data)<-pca_data[,idloc]
            #remove sampleids
            pca_data<-pca_data[-idloc]
            pca_data[pca_data==0]<-NA
            #remove empty columns
            empty_columns<-colSums(is.na(pca_data) | pca_data == "") == nrow(pca_data)
            pca_data<-pca_data[, !empty_columns]
            #replace NA with min per column
            pca_data<-pca_data %>%mutate_if(is.numeric,function(x) ifelse(is.na(x),min(x,na.rm=T),x))
            #remove constants
            pca_data<-remove_constant(pca_data)
            #pca
            nums<-unlist(lapply(pca_data,is.numeric))
            pca_obj<-prcomp(pca_data[,nums],center=TRUE,scale.=TRUE)
          })
          
          pca_table<-reactive({
            req(input$pca_file)
            pca<-fortify(pca_calc(),data=pca_file())
            dfcols<-ncol(pca_file())
            pca<-pca[,-(1:dfcols),drop=FALSE]
            pca = pca %>% `rownames<-`( NULL )
          })
          
          
          #pca_score
          output$pca_scores<-renderPlotly({
            req(input$pca_file,input$pca_id,input$pca_group)
            if(pca_labels()=="Yes"){
              pca_scores<-autoplot(pca_calc(), data=pca_file(),colour=pca_groups(),frame=TRUE,frame.type="norm",size=2,label=TRUE)+theme_classic()+theme(legend.position="bottom")
            } else {
              pca_scores<-autoplot(pca_calc(), data=pca_file(),colour=pca_groups(),frame=TRUE,frame.type="norm",size=2)+theme_classic()+theme(legend.position="bottom")
            }
            ggplotly(pca_scores) %>%
              layout(showlegend = TRUE, legend = list(font = list(size = 15)))
          })
          
          
          #pca_biplot
          output$pca_biplot<-renderPlotly({
            req(input$pca_file,input$pca_id,input$pca_group)
            pca_biplot<-autoplot(pca_calc(), data=pca_file(),colour=pca_groups(),loadings=TRUE,loadings.label=TRUE,loadings.label.size=3)+theme_classic()+theme(legend.position="none")
            ggplotly(pca_biplot) %>%
              layout(showlegend = TRUE, legend = list(font = list(size = 15)))
          })
          output$pca_table<-DT::renderDataTable({
            req(input$pca_file)
            pca_file()
          })
          output$pca_pca<-DT::renderDataTable({
            req(input$pca_file)
            pca_table()
            #row.names(pca)<-NULL
          })
          output$pca_download <- downloadHandler(
            filename = function() {
              paste("PCA_", ".csv", sep = "")
            },
            content = function(file) {
              write.csv(pca_table(), file, row.names = FALSE)
            }
          )
          
}

# Run the app ----
#shinyApp(ui = app_ui, server = app_server)
