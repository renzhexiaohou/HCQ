#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(ggplot2)
library(mrgsolve)
library(deSolve)
library(PKNCA)

code_pk <- '
        $PARAM KA = 1.15, CL = 10.9, VC = 437, VP = 1390, Q = 45.1

        $MAIN
        double CLi = CL*exp(ETA_CL);
        double VCi = VC*exp(ETA_VC);
        double VPi = VP*exp(ETA_VP);
        double Qi = Q;

        $INIT GUT = 0, CENT = 0, PERH = 0
        
        $ODE 
        double PKCONC = CENT/VCi*1000 / 433.9*336;
        dxdt_GUT = -KA*GUT;
        dxdt_CENT = KA*GUT - CLi/VCi*CENT - Qi/VCi*CENT + Qi/VPi*PERH;
        dxdt_PERH =  Qi/VCi*CENT - Qi/VPi*PERH;

        $OMEGA @name IIV @labels ETA_CL ETA_VC ETA_VP
        0 0 0
        
        $SIGMA 0
        
        $CAPTURE @annotated
        PKCONC: Concentration (ng/mL)
    '
mod_pk <- mcode("pk_model", code_pk)


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    
    dat_pk <- reactive({
        
        cyc <- input$cycle_pk
        
        amt1 <- input$amt_pk1 * input$p_pk
        tau1 <- input$interval_pk1
        n1 <- input$n_pk1
        
        amt2 <- input$amt_pk2 * input$p_pk
        tau2 <- input$interval_pk2
        n2 <- input$n_pk2
        
        cl <-input$cl_pk
        vc <- input$vc_pk
        vp <- input$vp_pk
        q <- input$q_pk
        ka <- input$ka_pk

        dos1 <- expand.ev(cmt = 1,
                          time = c(0,
                                   c(1:cyc) * (tau1*n1 + tau2*n2)),
                          amt = amt1,
                          ii = tau1,
                          addl = n1 - 1,
                          CL = cl,
                          VC = vc, 
                          VP = vp, 
                          Q = q, 
                          KA = ka,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_VC = rnorm(1, 0, 0),
                          ETA_VP = rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos2 <- expand.ev(cmt = 1,
                          time = c(c(tau1*n1 + tau2*n2, 
                                     c(1:(cyc)) * (tau1*n1 + tau2*n2)) - tau2*n2)[-1],
                          amt = amt2,
                          ii = tau2,
                          addl = n2 - 1,
                          CL = cl,
                          VC = vc, 
                          VP = vp, 
                          Q = q, 
                          KA = ka,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_VC = rnorm(1, 0, 0),
                          ETA_VP = rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos <- expand.ev(cmt = 1,
                         time = 0,
                         amt = 0,
                         ii = 24,
                         addl = 100,
                         CL = cl,
                         VC = vc, 
                         VP = vp, 
                         Q = q, 
                         KA = ka,
                         ETA_CL = rnorm(1, 0, 0), 
                         ETA_VC = rnorm(1, 0, 0),
                         ETA_VP = rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 2)
        
        # union_all(dos1, dos2) %>% arrange(ID, time, addl, ii)
        union_all(dos1, dos2) %>%
            # union_all(dos) %>% 
            arrange(ID, time, addl, ii)
    })
    
    output$distPlotPK <- renderPlot({
        cyc <- input$cycle_pk
        
        amt1 <- input$amt_pk1 * input$p_pk
        tau1 <- input$interval_pk1
        n1 <- input$n_pk1
        
        amt2 <- input$amt_pk2 * input$p_pk
        tau2 <- input$interval_pk2
        n2 <- input$n_pk2
        
        out <- mod_pk %>%
          data_set(dat_pk()) %>% 
          mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                 delta = 0.05)
        
        mod_pk %>%
            data_set(dat_pk()) %>% 
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) %>% 
            plot(PKCONC~time)
        
        pk <- out@data %>%
          select(time, PKCONC) %>% 
          mutate(subject = 1, time = time, conc = PKCONC) %>%
          distinct() %>%
          as.data.frame()
        
        ggplot(pk) +
          geom_line(aes(x = time/24, y = conc)) +
          geom_hline(yintercept = 21.1, linetype = "dashed", colour = "red") +
          scale_y_continuous(breaks = seq(0, 2000, 100)) +
          scale_x_continuous(breaks = seq(0, 7*100, 1)) +
          # scale_y_log10() +
          xlab("Time (Days)") + ylab("HCQ plasma concentration (ng/mL)")+
          theme_bw(base_size = 15)
        
        
    })
    
    output$halflife <- renderText({
        # browser()
        round(0.693/(input$cl_pk/input$v_pk), digits = 1)
    })
    
    
    reaccmax <- reactive({
        # PKNCA packages code as followed
        # (since ncappc packages can't be used in linux sys)
        
        cyc <- input$cycle_pk
        
        amt1 <- input$amt_pk1 * input$p_pk
        tau1 <- input$interval_pk1
        n1 <- input$n_pk1
        
        amt2 <- input$amt_pk2 * input$p_pk
        tau2 <- input$interval_pk2
        n2 <- input$n_pk2
        
        out <- mod_pk %>%
            data_set(dat_pk()) %>% 
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05)
        
        pk <- out@data %>%
            select(time, PKCONC) %>% 
            mutate(subject = 1, time = time, conc = PKCONC) %>%
            distinct() %>%
            as.data.frame()
        
        # browser()
        
        
        nca <- PKNCAconc(pk, conc~time|subject) %>%
            PKNCAdata(.,
                      intervals=data.frame(
                          start= cyc*(tau1*(n1) + tau2*(n2)) - tau1 - tau2*n2,
                          end= cyc*(tau1*(n1) + tau2*(n2)) - tau2*n2,
                          cmax=TRUE,
                          tmax=TRUE,
                          # aucinf.obs=TRUE,
                          auclast=TRUE
                      )
            ) %>%
            pk.nca() %>%
            as.data.frame(.$result)
        
        nca %>% filter(PPTESTCD == "cmax") %>% .$PPORRES %>% 
                round(digits = 1)
    })
    
    output$pkcmax1 <- renderText({
        # browser()
        reaccmax()
    })
    
    reacauc <- reactive({
        # PKNCA packages code as followed
        # (since ncappc packages can't be used in linux sys)
        
        cyc <- input$cycle_pk
        
        amt1 <- input$amt_pk1
        tau1 <- input$interval_pk1
        n1 <- input$n_pk1 * input$p_pk
        
        amt2 <- input$amt_pk2 * input$p_pk
        tau2 <- input$interval_pk2
        n2 <- input$n_pk2
        
        out <- mod_pk %>%
            data_set(dat_pk()) %>% 
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05)
        
        pk <- out@data %>%
            select(time, PKCONC) %>% 
            mutate(subject = 1, time = time, conc = PKCONC) %>%
            distinct() %>%
            as.data.frame()
        
        # browser()
        
        
        nca <- PKNCAconc(pk, conc~time|subject) %>%
            PKNCAdata(.,
                      intervals=data.frame(
                          start= cyc*(tau1*(n1) + tau2*(n2)) - tau1 - tau2*n2,
                          end= cyc*(tau1*(n1) + tau2*(n2)) - tau2*n2,
                          cmax=TRUE,
                          tmax=TRUE,
                          # aucinf.obs=TRUE,
                          auclast=TRUE
                      )
            ) %>%
            pk.nca() %>%
            as.data.frame(.$result)
        
        nca %>% filter(PPTESTCD == "auclast") %>% .$PPORRES %>% 
            round(digits = 1)
    })
    
    output$pkauc1 <- renderText({
        # browser()
        reacauc()
    })
    
  
})

