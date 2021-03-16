library(ggplot2)
library(reshape2)
library(scales)
#library(pracma)
library(doParallel)
library(deSolve)
library(cowplot)
library(magick)
library(tidyr)
library(viridis)

### PLOT OPTIONS ####

path_figure_results="Figures_results/"
path_figure_supp="Figures_supp_S2/"

theme<-theme_gray()+
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(colour='grey'),
        panel.grid.major.y = element_line(colour='grey'),
        text = element_text(size=20),
        axis.text = element_text(size=20),
        axis.line = element_line(),
        legend.key=element_blank(),
        plot.title = element_text(hjust = 0.5))

theme_matrix<-theme_gray()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text = element_text(size=20),
        axis.text = element_text(size=20),
        axis.line = element_line(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.key=element_blank(),
        plot.title = element_text(hjust = 0.5))

corr_colour_TL_2<-scale_colour_manual(values=c("dodgerblue3","chocolate1"),
                                      labels=c("1","2"),
                                      guide = guide_legend(reverse = TRUE),
                                      name='Trophic\nlevel')
corr_colour_TL_3<-scale_colour_manual(values=c("dodgerblue3","chocolate1","chartreuse4"),
                                      labels=c("1","2","3"),
                                      guide = guide_legend(reverse = TRUE),
                                      name='Trophic\nlevel')
corr_colour_TL_4<-scale_colour_manual(values=c("dodgerblue3","chocolate1","chartreuse4","red"),
                                      labels=c("1","2","3","4"),
                                      guide = guide_legend(reverse = TRUE),
                                      name='Trophic\nlevel')
corr_colour_TL_2_TS<-scale_colour_manual(values=c("dodgerblue1","chocolate1","dodgerblue4","chocolate4"),
                                         guide = guide_legend(reverse = TRUE),
                                         name='Trophic\nlevel')
model_line<-scale_linetype_manual(values=c("solid","twodash"),
                                  guide = guide_legend(reverse = TRUE),
                                  name='Dispersal')
disp_line<-scale_linetype_manual(values=c("solid","22"),
                                 labels=c("density\n-dependent","passive"),
                                 name='Dispersal')
patch_line<-scale_linetype_manual(values=c("solid","22"),
                                  name='Patch')
weight_line<-scale_linetype_manual(values=c("22","solid"),
                                  name='Weight of\ndispersal\ndependencies')
corr_line_TL_2_TS<-scale_linetype_manual(values=c("solid","solid","22","22"),
                                         guide = guide_legend(reverse = TRUE),
                                         name='Trophic\nlevel')
CV_colour_grad<-scale_fill_gradient2(low = "red",
                                     mid = "white",
                                     high = "blue",
                                     midpoint = 0,
                                     name="Covariance")
fill_colour_TL_4<-scale_fill_manual(values=c("dodgerblue3","chocolate1","chartreuse4","red"),
                                    labels=c("1","2","3","4"),
                                    guide = guide_legend(reverse = TRUE),
                                    name='Trophic\nlevel')
fill_colour_TL_3<-scale_fill_manual(values=c("dodgerblue3","chocolate1","chartreuse4"),
                                    labels=c("1","2","3"),
                                    guide = guide_legend(reverse = TRUE),
                                    name='Trophic\nlevel')
fill_colour_weight<-scale_fill_viridis(discrete=TRUE,
                                       name='Density\ndependency')
corr_colour_grad<-scale_fill_gradient2(low = "red",
                                       mid = "white",
                                       high = "blue",
                                       midpoint = 0,
                                       limits = c(-1,1),
                                       name="Correlation\ncoefficient")

x_axis_log10_short<-scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
y_axis_log10_short<-scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
x_axis_log10<-scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))
y_axis_log10<-scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))

# ea and ma factors
ma01<-expression(italic(ma)*"=0.1")
ma1<-expression(italic(ma)*"=1")
ma10<-expression(italic(ma)*"=10")
ea01<-expression(italic("\u03B5")*italic(a)*"=0.1")
ea1<-expression(italic("\u03B5")*italic(a)*"=1")
ea10<-expression(italic("\u03B5")*italic(a)*"=10")

# labels
label_dispersal<-expression("Scaled dispersal rate "*italic(d["i"]))
label_correlation<-"Correlation between the two patches"
label_S0<-expression(paste("Sensitivity coefficient ",italic(S["0"])))
label_CV<-"Coefficient of variation (CV)"
label_balance<-"Relative importance of dispersal"

### FUNCTIONS ####
# Define S values depending on B
set_S<-function(B,S0,nSpecies,nCommunity,type){
  S<-rep(NA,nSpecies*nCommunity)
  if(type=="self"){ # dispersal depending on the focal species
    for(i in 1:nSpecies){
      for(j in 1:nCommunity){
        S[(j-1)*nSpecies+i]=S0*B[(j-1)*nSpecies+i]
      }
    }
  }
  if(type=="prey"){ # dispersal depending on prey
    for(i in 2:nSpecies){
      for(j in 1:nCommunity){
        S[(j-1)*nSpecies+i]=S0*B[(j-1)*nSpecies+i-1]
      }
    }
  }
  if(type=="pred"){ # dispersal depending on predator
    for(i in 1:(nSpecies-1)){
      for(j in 1:nCommunity){
        S[(j-1)*nSpecies+i]=S0*B[(j-1)*nSpecies+i+1]
      }
    }
  }
  return(S)
}

# Correction of the dispersal rate to avoid the bias induced by density-dependency
set_C<-function(B,S,nSpecies,nCommunity,type,params){
  with(as.list(params),{
    # 0 - simple dispersal
    W<-rep(C0,nSpecies*nCommunity) # rescaling of the final dispersal rate to 1
    if(type=="0"){
      correct<-rep(1,nSpecies*nCommunity) # rescaling of the dispersal component to 1
      C<-rep(C0,nSpecies*nCommunity) # weight of the component
    }
    # self
    W<-W+C0_self*B # rescaling of the final dispersal rate to 1
    if(type=="self"){
      correct<-B/(B+S) # rescaling of the dispersal component to 1
      if(weight==1){
        C<-C0_self*B # weight of the component
      }
      else{C<-rep(C0_self,nSpecies*nCommunity)}
    }
    # prey
    Bbis<-B
    Bbis[2:(nSpecies*nCommunity)]<-Bbis[1:(nSpecies*nCommunity-1)]
    Bbis[is.na(S)]=0
    W<-W+C0_prey*e*a*Bbis # rescaling of the final dispersal rate to 1
    if(type=="prey"){
      correct<-S/(Bbis+S) # rescaling of the dispersal component to 1
      if(weight==1){
        C<-C0_prey*e*a*Bbis # weight of the component
      }
      else{C<-rep(C0_prey,nSpecies*nCommunity)}
    }
    # pred
    Bbis<-B
    Bbis[1:(nSpecies*nCommunity-1)]<-Bbis[2:(nSpecies*nCommunity)]
    Bbis[is.na(S)]=0
    W<-W+C0_pred*m*a*Bbis # rescaling of the final dispersal rate to 1
    if(type=="pred"){
      correct<-Bbis/(Bbis+S) # rescaling of the dispersal component to 1
      if(weight==1){
        C<-C0_pred*m*a*Bbis # weight of the component
      }
      else{C<-rep(C0_pred,nSpecies*nCommunity)}
    }
    # Final
    if(weight==1){
      C<-C/W # weighted components
    }
    # final correction
    correct[is.na(correct)]=1
    if(correction==1){ # rescaling of the dispersal component
      # effect of the different dispersal components
      if(nSpecies>1){
        prey<-rep(1,nSpecies*nCommunity)
        pred<-rep(1,nSpecies*nCommunity)
        for(j in 1:nCommunity){
          prey[(j-1)*nSpecies+1]=0 # no prey for primary producers
          pred[(j-1)*nSpecies+nSpecies]=0 # ni predator for top predators
        }
      }
      else{prey<-rep(0,nSpecies*nCommunity)
      pred<-rep(0,nSpecies*nCommunity)}
      if(weight==0){ # dispersal component not weighted by the intra-patch processes
        C<-C/(C0+C0_self+prey*C0_prey+pred*C0_pred) # the dispersal components have the same weight (because C0=1, C0_self=1...)
        C[is.na(C)]=0
      }
      C<-C/correct # correction to get a final dispersal rate with a similar value than the case with passive dispersal
    }
    return(C)
  })
}

# Jacobian single patch J(k)
jacobian_single<-function(B,params,nSpecies){
  with(as.list(params),{
    J<-matrix(0,nrow=nSpecies,ncol=nSpecies)
    J[1,1]=D*(g/D-2*B[1])
    if(nSpecies>1){
      J[1,1]=J[1,1]-D*m*a*B[2] # basal species
      J[nSpecies,nSpecies]=D*m^(nSpecies-1)*(-r/D-2*B[nSpecies]+e*a*B[nSpecies-1]) # top species
      for(i in 1:(nSpecies-1)){
        J[i,i+1]=-D*m^(i-1)*m*a*B[i] # effect of predators
      }
      for(i in 2:nSpecies){
        J[i,i-1]=D*m^(i-1)*e*a*B[i] # effect of prey
      }
    }
    if(nSpecies>2){
      for(i in 2:(nSpecies-1)){
        J[i,i]=D*m^(i-1)*(-r/D-2*B[i]+e*a*B[i-1]-m*a*B[i+1]) # intermediate species
      }
    }
    return(J)
  })
}

# Jacobian of the intrapatch dynamics J'
jacobian_intra_patch<-function(B,params,nSpecies,nCommunity){
  Z<-matrix(0,nrow=nSpecies,ncol=nSpecies)
  M<-NULL
  L<-NULL
  J<-NULL
  for(i in 1:nCommunity){
    L<-NULL
    J<-jacobian_single(B[nSpecies*(i-1)+c(1:nSpecies)],params,nSpecies)
    for(j in 1:nCommunity){
      if(j==i){
        L<-cbind(L,J)
      }
      else{
        L<-cbind(L,Z)
      }
    }
    M<-rbind(M,L)
  }
  return(M)
}

# Jacobian of the dispersal dynamics P'
jacobian_disp_hauzy<-function(B,params,nSpecies,nCommunity,disp){
  with(as.list(params),{
    dim=nSpecies*nCommunity
    P<-matrix(0,nrow=dim,ncol=dim)
    # effects of species i on itself
    S=set_S(B,S0_self,nSpecies,nCommunity,"self")
    C_0=set_C(B,NA,nSpecies,nCommunity,"0",params)
    C_self=set_C(B,S,nSpecies,nCommunity,"self",params)
    for(i in 1:nSpecies){
      for(j in 1:nCommunity){
        for(k in 1:nCommunity){
          row=(j-1)*nSpecies+i
          col=(k-1)*nSpecies+i
          if(j==k){
            P[row,col]=-disp[i]*(C_0[col]+C_self[col]*B[col]*(B[col]+2*S[col])/(B[col]+S[col])^2)
          }
          else{
            P[row,col]=disp[i]*(C_0[col]+C_self[col]*B[col]*(B[col]+2*S[col])/(B[col]+S[col])^2)/(nCommunity-1)
          }
        }
      }
    }
    # effects of prey
    S=set_S(B,S0_prey,nSpecies,nCommunity,"prey")
    C=set_C(B,S,nSpecies,nCommunity,"prey",params)
    for(i in 2:nSpecies){
      for(j in 1:nCommunity){
        for(k in 1:nCommunity){
          row=(j-1)*nSpecies+i
          col=(k-1)*nSpecies+i
          if(j==k){
            P[row,col]=P[row,col]-disp[i]*C[col]*S[col]/(B[col-1]+S[col])
            P[row,col-1]=disp[i]*B[col]*C[col]*S[col]/(B[col-1]+S[col])^2
          }
          else{
            P[row,col]=P[row,col]+disp[i]*C[col]*S[col]/(B[col-1]+S[col])/(nCommunity-1)
            P[row,col-1]=-disp[i]*B[col]*C[col]*S[col]/(B[col-1]+S[col])^2/(nCommunity-1)
          }
        }
      }
    }
    # effects of predators
    S=set_S(B,S0_pred,nSpecies,nCommunity,"pred")
    C=set_C(B,S,nSpecies,nCommunity,"pred",params)
    for(i in 1:(nSpecies-1)){
      for(j in 1:nCommunity){
        for(k in 1:nCommunity){
          row=(j-1)*nSpecies+i
          col=(k-1)*nSpecies+i
          if(j==k){
            P[row,col]=P[row,col]-disp[i]*C[col]*B[col+1]/(B[col+1]+S[col])
            P[row,col+1]=-disp[i]*B[col]*C[col]*S[col]/(B[col+1]+S[col])^2
          }
          else{
            P[row,col]=P[row,col]+disp[i]*C[col]*B[col+1]/(B[col+1]+S[col])/(nCommunity-1)
            P[row,col+1]=disp[i]*B[col]*C[col]*S[col]/(B[col+1]+S[col])^2/(nCommunity-1)
          }
        }
      }
    }
    # Finalising the matrix
    for(i in 2:nSpecies){
      for(j in 1:nCommunity){
        P[(j-1)*nSpecies+i,]=D*m^(i-1)*P[(j-1)*nSpecies+i,]
      }
    }
    return(P)
  })
}

# Jacobian of the dispersal dynamics P' (linear dispersal components)
jacobian_disp_linear<-function(B,params,nSpecies,nCommunity,disp){
  with(as.list(params),{
    dim=nSpecies*nCommunity
    P<-matrix(0,nrow=dim,ncol=dim)
    correct<-rep(0,dim)
    # effects of species i on itself
    if(correction==1){ # rescaling of the dispersal component
      C<-C0_self/B
      correct<-correct+C0_self
    }else{C<-rep(C0_self,dim)}
    for(i in 1:nSpecies){
      for(j in 1:nCommunity){
        for(k in 1:nCommunity){
          row=(j-1)*nSpecies+i
          col=(k-1)*nSpecies+i
          if(j==k){
            P[row,col]=-disp[i]*2*C[col]*B[col]
          }
          else{
            P[row,col]=disp[i]*2*C[col]*B[col]
          }
        }
      }
    }
    # effects of prey
    if(correction==1){ # rescaling of the dispersal component
      C<-C0_prey/c(0,B[1:(dim-1)])
      C[seq(from=1,to=(nSpecies*(nCommunity-1)+1),by=nSpecies)]=0
      correct<-correct+C+C0
    }else{C<-rep(C0_prey,dim)}
    for(i in 2:nSpecies){
      for(j in 1:nCommunity){
        for(k in 1:nCommunity){
          row=(j-1)*nSpecies+i
          col=(k-1)*nSpecies+i
          if(j==k){
            P[row,col]=P[row,col]-disp[i]*(-C[col]*B[col-1]+C0)
            P[row,col-1]=disp[i]*C[col]*B[col]
          }
          else{
            P[row,col]=P[row,col]+disp[i]*(-C[col]*B[col-1]+C0)
            P[row,col-1]=-disp[i]*C[col]*B[col]
          }
        }
      }
    }
    # effects of predators
    if(correction==1){ # rescaling of the dispersal component
      C<-C0_pred/c(B[2:dim],NA)
      C[seq(from=nSpecies,to=(nSpecies*nCommunity),by=nSpecies)]=0
      correct<-correct+C
    }else{C<-rep(C0_pred,dim)}
    for(i in 1:(nSpecies-1)){
      for(j in 1:nCommunity){
        for(k in 1:nCommunity){
          row=(j-1)*nSpecies+i
          col=(k-1)*nSpecies+i
          if(j==k){
            P[row,col]=P[row,col]-disp[i]*(C[col]*B[col+1])
            P[row,col+1]=-disp[i]*C[col]*B[col]
          }
          else{
            P[row,col]=P[row,col]+disp[i]*(C[col]*B[col+1])
            P[row,col+1]=disp[i]*C[col]*B[col]
          }
        }
      }
    }
    # Finalising the matrix
    for(i in 2:nSpecies){
      for(j in 1:nCommunity){
        P[(j-1)*nSpecies+i,]=D*m^(i-1)*P[(j-1)*nSpecies+i,]
      }
    }
    if(correction==1){ # rescaling of the dispersal component
      # effect of the different dispersal components
      if(nSpecies>1){
        prey<-rep(1,nSpecies*nCommunity)
        pred<-rep(1,nSpecies*nCommunity)
        for(j in 1:nCommunity){
          prey[(j-1)*nSpecies+1]=0 # no prey for primary producers
          pred[(j-1)*nSpecies+nSpecies]=0 # ni predator for top predators
        }
      }
      else{prey<-rep(0,nSpecies*nCommunity)
      pred<-rep(0,nSpecies*nCommunity)}
      if(weight==0){ # dispersal component not weighted by the intra-patch processes
        for(i in 1:dim){
          P[i,]=P[i,]/(C0+C0_self+prey*C0_prey+pred*C0_pred) # final correction to rescale the dispersal rate
        }
        P[is.na(P)]=0
      }
    }
    return(P)
  })
}

# Biomasses at equilibrium
equilibrium<-function(params,nSpecies,nCommunity){ # compute the biomasses at equilibrium
  with(as.list(params),{
    A<-diag(rep(-1,nSpecies*nCommunity))
    for(j in 1:nCommunity){
      for(i in 2:nSpecies){
        A[(j-1)*nSpecies+i,(j-1)*nSpecies+i-1]=e*a
      }
      for(i in 1:(nSpecies-1)){
        A[(j-1)*nSpecies+i,(j-1)*nSpecies+i+1]=-m*a
      }
    }
    B<-matrix(r/D,
              nrow = nSpecies*nCommunity,
              ncol = 1,
              byrow = TRUE)
    for(j in 1:nCommunity){
      B[(j-1)*nSpecies+1,1]=-g/D
    }
    C<-matrix(0,
              nrow = nSpecies*nCommunity,
              ncol = 1,
              byrow = TRUE)
    C<-solve(A) %*% B
    return(as.numeric(C))
  })
}

# T matrix
T_matrix<-function(pert,B,z,nSpecies,nCommunity){
  T<-matrix(0,nSpecies*nCommunity,nSpecies*nCommunity)
  coord=NULL
  for(i in 1:length(pert)){
    coord=(pert[[i]][2]-1)*nSpecies+pert[[i]][1] # pert->list of vector containing the trophic level and the patch of the perturbed species (e.g. (1,2) species 1 in patch 2)
    T[coord,coord]=1
  }
  T<-T*diag(B)^z
  return(T)
}

# Lyapunov equation
lyapunov<-function(J,T,VE,nSpecies,nCommunity){
  TVT<-T%*%VE%*%t(T)
  TVT<-matrix(array(TVT),ncol=1)
  kron<-kronecker(J,diag(rep(1,nSpecies*nCommunity))) + kronecker(diag(rep(1,nSpecies*nCommunity)),J)
  return(-solve(kron)%*%TVT)
}

# Parallelised analytic calculation 
analytic<-function(params_data,jacobian_disp,nSpecies,nCommunity,i){
  params<-c(g=params_data$g[i],
            r=params_data$r[i],
            D=params_data$D[i],
            m=params_data$m[i],
            a=params_data$a[i],
            e=params_data$e[i],
            S0_self=params_data$S0_self[i],
            S0_prey=params_data$S0_prey[i],
            S0_pred=params_data$S0_pred[i],
            C0=params_data$C0[i],
            C0_self=params_data$C0_self[i],
            C0_prey=params_data$C0_prey[i],
            C0_pred=params_data$C0_pred[i],
            correction=params_data$correction[i],
            weight=params_data$weight[i])
  B<-equilibrium(params,nSpecies,nCommunity)
  J<-jacobian_intra_patch(B,params,nSpecies,nCommunity)
  P<-jacobian_disp(B,params,nSpecies,nCommunity,params_data$disp[[i]]) # jacobian_disp_hauzy or jacobian_disp_linear
  VE<-params_data$VE[[i]]
  T<-T_matrix(params_data$pert[[i]],B,params_data$z[i],nSpecies,nCommunity)
  V<-lyapunov(J+P,T,VE,nSpecies,nCommunity)
  C<-cov2cor(matrix(as.numeric(V),nSpecies*nCommunity,nSpecies*nCommunity))
  return(list(B,as.numeric(V),as.numeric(C)))
}

# Create the output dataframe
create_data<-function(params_data,nSpecies,nCommunity){
  # biomass
  data_B<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+nSpecies*nCommunity))
  names(data_B)[1:dim(params_data)[2]]=names(params_data)
  data_B[,1:dim(params_data)[2]]=params_data
  B_names<-expand.grid(c("B"),c(1:nSpecies),c(1:nCommunity))
  B_names$Var1<-paste(B_names$Var1,B_names$Var2,sep = "_")
  B_names<-paste(B_names$Var1,B_names$Var3,sep = "")
  names(data_B)[(dim(params_data)[2]+1):dim(data_B)[2]]=B_names
  # variance
  data_V<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+(nSpecies*nCommunity)^2))
  names(data_V)[1:dim(params_data)[2]]=names(params_data)
  data_V[,1:dim(params_data)[2]]=params_data
  V_names<-expand.grid(c(1:nSpecies),c(1:nCommunity))
  V_names<-paste(V_names$Var1,V_names$Var2,sep = "")
  V_names<-expand.grid(c("V"),V_names,V_names)
  V_names<-paste(V_names$Var1,V_names$Var2,V_names$Var3,sep = "_")
  names(data_V)[(dim(params_data)[2]+1):dim(data_V)[2]]=V_names
  # correlation
  data_C<-as.data.frame(matrix(0,dim(params_data)[1],dim(params_data)[2]+(nSpecies*nCommunity)^2))
  names(data_C)[1:dim(params_data)[2]]=names(params_data)
  data_C[,1:dim(params_data)[2]]=params_data
  C_names<-expand.grid(c(1:nSpecies),c(1:nCommunity))
  C_names<-paste(C_names$Var1,C_names$Var2,sep = "")
  C_names<-expand.grid(c("C"),C_names,C_names)
  C_names<-paste(C_names$Var1,C_names$Var2,C_names$Var3,sep = "_")
  names(data_C)[(dim(params_data)[2]+1):dim(data_C)[2]]=C_names
  return(list(data_B,data_V,data_C,B_names,V_names,C_names))
}

# ODE of the system
ODE_function<-function(t, B, params){
  with(params,{
    dB<-rep(0,nSpecies*nCommunity)
    # intra-patch dynamics
    ID=0
    for(j in 1:nCommunity){
      ID=(j-1)*nSpecies
      dB[ID+1]=D*B[ID+1]*(g/D - B[ID+1])
      if(nSpecies>1){
        for(i in 2:nSpecies){
          dB[ID+i]=D*m^(i-1)*B[ID+i]*(-r/D - B[ID+i] + e*a*B[ID+i-1])
        }
        for(i in 1:(nSpecies-1)){
          dB[ID+i]=dB[ID+i]+D*m^(i-1)*B[ID+i]*(- m*a*B[ID+i+1])
        }
      }
    }
    # dispersal dynamics
    d_matrix<-matrix(C0,nSpecies,nCommunity) # matrix containing the dispersal terms
    B_matrix<-matrix(B,nSpecies,nCommunity) # matrix containg the biomasses
    # effects of species i on itself
    S_matrix<-matrix(S_self,nSpecies,nCommunity) # matrix containing the dispersal sensitivity S
    C_matrix<-matrix(C_self,nSpecies,nCommunity) # matrix containing the correction coefficients
    for(i in 1:nSpecies){
      d_matrix[i,]=d_matrix[i,]+C_matrix[i,]*B_matrix[i,]/(B_matrix[i,]+S_matrix[i,])
    }
    # effects of prey
    S_matrix<-matrix(S_prey,nSpecies,nCommunity)
    C_matrix<-matrix(C_prey,nSpecies,nCommunity)
    if(nSpecies>1){
      for(i in 2:nSpecies){
        d_matrix[i,]=d_matrix[i,]+C_matrix[i,]*S_matrix[i,]/(B_matrix[i-1,]+S_matrix[i,])
      }
    }
    # effects of predators
    S_matrix<-matrix(S_pred,nSpecies,nCommunity)
    C_matrix<-matrix(C_pred,nSpecies,nCommunity)
    if(nSpecies>1){
      for(i in 1:(nSpecies-1)){
        d_matrix[i,]=d_matrix[i,]+C_matrix[i,]*B_matrix[i+1,]/(B_matrix[i+1,]+S_matrix[i,])
      }
    }
    # final equations
    for(i in 1:nSpecies){
      d_matrix[i,]=m^(i-1)*D*disp[i]*B_matrix[i,]*d_matrix[i,]
    }
    for(i in 1:nSpecies){
      for(j in 1:nCommunity){
        dB[(j-1)*nSpecies+i]=dB[(j-1)*nSpecies+i]-(1+1/(nCommunity-1))*d_matrix[i,j]+sum(d_matrix[i,])/(nCommunity-1)
      }
    }
    return(list(dB))
  })
}

# Parallelised function solving the ODE
ODE_solve<-function(time, params_data, nSpecies, nCommunity, i){
  params<-c(g=params_data$g[i],
            r=params_data$r[i],
            D=params_data$D[i],
            m=params_data$m[i],
            a=params_data$a[i],
            e=params_data$e[i],
            nSpecies=nSpecies,
            nCommunity=nCommunity)
  B0<-equilibrium(params,nSpecies,nCommunity)
  S_self=set_S(B0,params_data$S0_self[i],nSpecies,nCommunity,"self")
  S_prey=set_S(B0,params_data$S0_prey[i],nSpecies,nCommunity,"prey")
  S_pred=set_S(B0,params_data$S0_prey[i],nSpecies,nCommunity,"pred")
  C_0=set_C(B0,NA,nSpecies,nCommunity,"0",params_data)
  C_self=set_C(B0,S_self,nSpecies,nCommunity,"self",params_data)
  C_prey=set_C(B0,S_prey,nSpecies,nCommunity,"prey",params_data)
  C_pred=set_C(B0,S_pred,nSpecies,nCommunity,"pred",params_data)
  params<-c(as.list(params),
            list(S_self=S_self),list(S_prey=S_prey),list(S_pred=S_pred),
            list(C_self=C_self),list(C_prey=C_prey),list(C_pred=C_pred),
            list(disp=params_data$disp[[i]]))
  for(j in length(params_data$pert[[i]])){
    sp<-params_data$pert[[i]][[j]]
    B0[(sp[2]-1)*nSpecies+sp[1]]=5*B0[(sp[2]-1)*nSpecies+sp[1]] # perturbation of a species
  }
  #B0[params_data$pert[[i]][1]]=5*B0[params_data$pert[[i]][1]]
  TS<-as.data.frame(ode(B0, time, ODE_function, params, method="rk4"))
  TS<-rbind(c(0,equilibrium(params,nSpecies,nCommunity)),TS) # add the biomasses at equilibrium to the time series
  return(TS)
}

# Make a table to plot a correlation matrix
table_for_matrix<-function(table,nparams){
  table<-melt(table,
              id.vars = names(table)[1:nparams],
              variable.name = "species",
              value.name = "value")
  table<-table %>% separate(species,c(NA,"species_1","species_2"),sep="_")
  table$species_1<-as.factor(table$species_1)
  table$species_2<-as.factor(table$species_2)
  table<-table %>% separate(species_1,c("species_1","community_1"),sep=1)
  table<-table %>% separate(species_2,c("species_2","community_2"),sep=1)
  return(table)
}

# Make a table to plot the biomass
table_for_biomass<-function(table,nparams){
  table<-melt(table,
              id.vars = names(table)[1:nparams],
              variable.name = "species",
              value.name = "biomass")
  table<-table %>% separate(species,int=c(NA,"species"),by="_")
  table$species<-as.factor(table$species)
  table<-table %>% separate(species,c("species","community"),sep=1)
  return(table)
}

# Retun the weight of each dispersal dependency
get_weight<-function(B,e,a,m,C0_self,C0_prey,C0_pred,nSpecies){
  W<-matrix(0,nrow=3,ncol=nSpecies)
  W[1,]<-C0_self*B
  if(nSpecies>1){
    W[2,2:nSpecies]<-C0_prey*e*a*B[1:(nSpecies-1)]
  }
  if(nSpecies>1){
    W[3,1:(nSpecies-1)]<-C0_pred*m*a*B[2:nSpecies]
  }
  for(i in 1:nSpecies){
    W[,i]<-W[,i]/sum(W[,i])
  }
  W<-as.data.frame(W)
  names(W)<-c(1:nSpecies)
  W$dependency<-c("self","prey","pred")
  W<-melt(W, id.vars = c("dependency"),
          variable.name = "species",
          value.name = "weight")
  return(W)
}

# dispersal importance for two patches
get_dispersal_importance<-function(params_data,nSpecies,i){
  params<-c(g=params_data$g[i],
            r=params_data$r[i],
            D=params_data$D[i],
            m=params_data$m[i],
            a=params_data$a[i],
            e=params_data$e[i],
            S0_self=params_data$S0_self[i],
            S0_prey=params_data$S0_prey[i],
            S0_pred=params_data$S0_pred[i],
            C0=params_data$C0[i],
            C0_self=params_data$C0_self[i],
            C0_prey=params_data$C0_prey[i],
            C0_pred=params_data$C0_pred[i],
            correction=params_data$correction[i],
            weight=params_data$weight[i],
            d=params_data$d[i])
  B<-equilibrium(params,nSpecies,1)
  # effects of trophic interactions
  trophic<-with(as.list(params),{
    trophic<-B
    trophic[1]=trophic[1]+g/D
    if(nSpecies>1){
      trophic[2:nSpecies]=trophic[2:nSpecies]+r/D+e*a*B[1:(nSpecies-1)]
      trophic[1:(nSpecies-1)]=trophic[1:(nSpecies-1)]+m*a*B[2:nSpecies]
    }
    trophic<-trophic*B
  })
  # effect of dispersal
  disp<-with(as.list(params),{
    # effects of basal dispersal
    disp<-set_C(B,NA,nSpecies,1,"0",params)
    # effects of self sensity dependency
    S=set_S(B,S0_self,nSpecies,1,"self")
    C=set_C(B,S,nSpecies,1,"self",params)
    disp<-disp+C*B/(B+S)
    if(nSpecies>1){
    # effects of prey
    S=set_S(B,S0_prey,nSpecies,1,"prey")
    C=set_C(B,S,nSpecies,1,"prey",params)
    C[is.na(S)]=0
    S[is.na(S)]=1
    Bbis<-B
    Bbis[2:nSpecies]<-B[1:(nSpecies-1)]
    disp<-disp+C*S/(Bbis+S)
    # effects of predators
    S=set_S(B,S0_pred,nSpecies,1,"pred")
    C=set_C(B,S,nSpecies,1,"pred",params)
    C[is.na(S)]=0
    S[is.na(S)]=1
    Bbis<-B
    Bbis[1:(nSpecies-1)]<-B[2:nSpecies]
    disp<-disp+C*Bbis/(Bbis+S)
    }
    disp<-2*d*B*disp
  })
  return(disp/(trophic+disp))
}

### PARAMETERS ####
d_min=-5
d_max=5
d_step=0.1
d=10^(seq(d_min,d_max,d_step))

S0_min=-4
S0_max=4
S0_step=0.1
S0=10^(seq(S0_min,S0_max,S0_step))

g=1
r=0
D=1
e=0.65
m=c(0.0065,0.065,0.65,6.5,65)
a=c(1/6.5,1/0.65,1/0.065)
sigma=1e-3
z=0.5
C0=0
C0_self=1
C0_prey=1
C0_pred=1
correction=1
weight=0
pert=list(c(1,1)) # c(species, patch)
disp=list(c(0,1))# c(species1, species2,..., nSpecies)

params_data_original<-expand.grid(simu_ID=0,
                                  g=g,
                                  r=r,
                                  D=D,
                                  e=e,
                                  m=m,
                                  a=a,
                                  sigma=sigma,
                                  z=z,
                                  weight=weight)
params_data_original$ma=params_data_original$m*params_data_original$a
params_data_original<-params_data_original[params_data_original$ma>0.05 & params_data_original$ma<=15,]
params_data_original$ea=params_data_original$e*params_data_original$a
params_data_original$simu_ID<-c(1:dim(params_data_original)[1])
params_data_original$ea<-as.factor(params_data_original$ea)
params_data_original$ma<-as.factor(params_data_original$ma)
levels(params_data_original$ma)<-c(ma01,ma1,ma10)
params_data_original$ma = factor(params_data_original$ma,levels(params_data_original$ma)[c(3,2,1)])
levels(params_data_original$ea)<-c(ea01,ea1,ea10)

########################## ----
# MAIN TEXT ############## ----
########################## ----
### RESULTS ############## ----
########################## ----
# DISPERSAL OF PRED DEPENDING ON PREY #### ----
### SIMULATIONS - d - CORRELATION ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         d=d,
                         S0_self=1,
                         S0_prey=c(1e-3,1e3),
                         S0_pred=1,
                         C0=0,
                         C0_self=0,
                         C0_prey=1,
                         C0_pred=0,
                         correction=correction,
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
params_data<-params_data[which(params_data$ea%in%as.character(c(ea01,ea10))),]
params_data<-params_data[which(params_data$ma%in%as.character(c(ma01,ma10))),]
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_prey","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","S0_prey"),
              variable.name = "Species",
              value.name = "Correlation")
databis$S0_prey<-as.factor(databis$S0_prey)

p1<-ggplot(data=databis)+
  geom_line(aes(S0_prey,Correlation,colour=Species,linetype=S0_prey),size=1.5)+
  corr_colour_TL_2+
  disp_line+
  theme
legend_1<-get_legend(p1)

p1<-ggplot(data=databis)+
      geom_vline(xintercept=1e3,linetype='dashed')+
      geom_line(aes(d,Correlation,colour=Species,linetype=S0_prey),size=1.5)+
      annotate("label",x=1e3,y=0.2,label="C",size=6)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      disp_line+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle(expression(atop("Dispersal of predators","depending on prey ("*italic(S["0"])*"=10"^"-3"*")")))

### SIMULATIONS - S0 - CORRELATION ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         d=1e3,
                         S0_self=1,
                         S0_prey=S0,
                         S0_pred=1,
                         C0=0,
                         C0_self=0,
                         C0_prey=1,
                         C0_pred=0,
                         correction=correction,
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
params_data<-params_data[which(params_data$ea%in%as.character(c(ea01,ea10))),]
params_data<-params_data[which(params_data$ma%in%as.character(c(ma01,ma10))),]
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_prey","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","S0_prey"),
              variable.name = "Species",
              value.name = "Correlation")

p3<-ggplot(data=databis[databis$d==1e3,])+
      geom_line(aes(S0_prey,Correlation,colour=Species),size=1.5)+
      corr_colour_TL_2+
      theme
legend_2<-get_legend(p3)

p3<-ggplot(data=databis[databis$d==1e3,])+
      geom_vline(xintercept=1e-3,linetype='dashed')+
      geom_line(aes(S0_prey,Correlation,colour=Species),size=1.5)+
      annotate("label",x=1e-3,y=0.2,label="A",size=6)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_S0)+
      ylab(label_correlation)+
      ggtitle(expression(atop("Dispersal of predators","depending on prey ("*italic(d["i"])*"=10"^"3"*")")))

# DISPERSAL OF PREY DEPENDING ON PRED #### ----
### SIMULATIONS - d - CORRELATION ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(2,1))), # perturbation of predators in patch 1
                         disp=list(c(1,0)), # dispersal of prey
                         d=d,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=c(1e-3,1e3),
                         C0=0,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=1,
                         correction=correction,
                         model="disp_prey")
params_data<-merge(params_data_original,params_data)
params_data<-params_data[which(params_data$ea%in%as.character(c(ea01,ea10))),]
params_data<-params_data[which(params_data$ma%in%as.character(c(ma01,ma10))),]
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_pred","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","S0_pred"),
              variable.name = "Species",
              value.name = "Correlation")
databis$S0_pred<-as.factor(databis$S0_pred)
databis$S0_pred = factor(databis$S0_pred,levels(databis$S0_pred)[c(2,1)])

p2<-ggplot(data=databis)+
      geom_vline(xintercept=1e3,linetype='dashed')+
      geom_line(aes(d,Correlation,colour=Species,linetype=S0_pred),size=1.5)+
      annotate("label",x=1e3,y=-0.5,label="D",size=6)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      disp_line+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle(expression(atop("Dispersal of prey","depending on predators ("*italic(S["0"])*"=10"^"3"*")")))

### SIMULATIONS - S0 - CORRELATION ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(2,1))), # perturbation of predators in patch 1
                         disp=list(c(1,0)), # dispersal of prey
                         d=1e3,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=S0,
                         C0=0,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=1,
                         correction=correction,
                         model="disp_prey")
params_data<-merge(params_data_original,params_data)
params_data<-params_data[which(params_data$ea%in%as.character(c(ea01,ea10))),]
params_data<-params_data[which(params_data$ma%in%as.character(c(ma01,ma10))),]
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_pred","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","S0_pred"),
              variable.name = "Species",
              value.name = "Correlation")

p4<-ggplot(data=databis[databis$d==1e3,])+
      geom_vline(xintercept=1e3,linetype='dashed')+
      geom_line(aes(S0_pred,Correlation,colour=Species),size=1.5)+
      annotate("label",x=1e3,y=-0.5,label="B",size=6)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_S0)+
      ylab(label_correlation)+
      ggtitle(expression(atop("Dispersal of prey","depending on predators ("*italic(d["i"])*"=10"^"3"*")")))

### FINAL GRAPH ####
p1_img<-image_read_pdf(paste(path_figure_results,"schema_little_density_pert_TL1_disp_TL2.pdf",sep=""))
p2_img<-image_read_pdf(paste(path_figure_results,"schema_little_density_pert_TL2_disp_TL1.pdf",sep=""))

graph<-ggdraw(xlim = c(0, 2.25), ylim = c(0, 2.35)) +
  draw_image(p1_img, 0.25, 2.05, 0.5, 0.3)+
  draw_image(p2_img, 1.25, 2.05, 0.5, 0.3)+
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(p4, 1, 0, 1, 1)+
  draw_plot(legend_1, 2.1, 1.25, 0.05, 0.5)+
  draw_plot(legend_2, 2.05, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 30)
ggsave(paste(path_figure_results,"figure_prey_pred.pdf",sep=""),graph, width = 15, height = 16, device = cairo_pdf)

########################## ----
# TOP-DOWN CONTROL #### ----
### SIMULATIONS - S0 - CORRELATION - BIOMASS ####
# parameters - 2 species
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(2,1))), # perturbation of predators in patch 1
                         disp=list(c(1,0)), # dispersal of prey
                         d=1e3,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=S0,
                         C0=0,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=1,
                         correction=correction,
                         model="2")
params_data<-merge(params_data_original,params_data)
params_data<-params_data[which(params_data$ea%in%as.character(c(ea01,ea10))),]
params_data<-params_data[which(params_data$ma%in%as.character(c(ma01,ma10))),]
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations - 2 species
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_pred","model","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","S0_pred","model"),
              variable.name = "Species",
              value.name = "Correlation")

data_B$B_31=NA
data_B$B_41=NA
databis2<-melt(data_B[,which(names(data_B)%in%c("ea","ma","d","S0_pred","model","B_11",'B_21',"B_31","B_41"))],
               id.vars = c("ea","ma","d","S0_pred","model"),
               variable.name = "Species",
               value.name = "Biomass")

# parameters - 3 species
nSpecies=3
nCommunity=2
params_data<-expand.grid(pert=list(list(c(2,1))), # perturbation of predators in patch 1
                         disp=list(c(1,0,0)), # dispersal of prey
                         d=1e3,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=S0,
                         C0=0,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=1,
                         correction=correction,
                         model="3")
params_data<-merge(params_data_original,params_data)
params_data<-params_data[which(params_data$ea%in%as.character(c(ea01,ea10))),]
params_data<-params_data[which(params_data$ma%in%as.character(c(ma01,ma10))),]
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations - 3 species
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-rbind(databis,melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_pred","model","C_11_12",'C_21_22'))],
                            id.vars = c("ea","ma","d","S0_pred","model"),
                            variable.name = "Species",
                            value.name = "Correlation"))

data_B$B_41=NA
databis2<-rbind(databis2,melt(data_B[,which(names(data_B)%in%c("ea","ma","d","S0_pred","model","B_11",'B_21',"B_31","B_41"))],
                              id.vars = c("ea","ma","d","S0_pred","model"),
                              variable.name = "Species",
                              value.name = "Biomass"))

# parameters - 4 species
nSpecies=4
nCommunity=2
params_data<-expand.grid(pert=list(list(c(2,1))), # perturbation of predators in patch 1
                         disp=list(c(1,0,0,0)), # dispersal of prey
                         d=1e3,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=S0,
                         C0=0,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=1,
                         correction=correction,
                         model="4")
params_data<-merge(params_data_original,params_data)
params_data<-params_data[which(params_data$ea%in%as.character(c(ea01,ea10))),]
params_data<-params_data[which(params_data$ma%in%as.character(c(ma01,ma10))),]
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations - 4 species
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-rbind(databis,databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_pred","model","C_11_12",'C_21_22'))],
                                     id.vars = c("ea","ma","d","S0_pred","model"),
                                     variable.name = "Species",
                                     value.name = "Correlation"))

databis2<-rbind(databis2,melt(data_B[,which(names(data_B)%in%c("ea","ma","d","S0_pred","model","B_11",'B_21',"B_31","B_41"))],
                              id.vars = c("ea","ma","d","S0_pred","model"),
                              variable.name = "Species",
                              value.name = "Biomass"))

databis2$x<-as.numeric(databis2$Species)
p1<-ggplot(data=databis2[databis2$ea==as.character(ea10) & databis2$ma==as.character(ma10) & databis2$model=="2",])+
      geom_rect(aes(xmin=x-0.5,xmax=x+0.5,ymin=1e-3,ymax=Biomass,fill=Species))+
      facet_grid(ma~ea, labeller=label_parsed)+
      fill_colour_TL_4+
      y_axis_log10_short+
      coord_flip()+
      theme+theme(legend.position = "none")+
      xlab('Trophic level')+
      ylab("Biomass")

p2<-ggplot(data=databis2[databis2$ea==as.character(ea10) & databis2$ma==as.character(ma10) & databis2$model=="3",])+
      geom_rect(aes(xmin=x-0.5,xmax=x+0.5,ymin=1e-3,ymax=Biomass,fill=Species))+
      facet_grid(ma~ea, labeller=label_parsed)+
      fill_colour_TL_4+
      y_axis_log10_short+
      coord_flip()+
      theme+theme(legend.position = "none")+
      xlab('Trophic level')+
      ylab("Biomass")

p3<-ggplot(data=databis2[databis2$ea==as.character(ea10) & databis2$ma==as.character(ma10) & databis2$model=="4",])+
      geom_rect(aes(xmin=x-0.5,xmax=x+0.5,ymin=1e-3,ymax=Biomass,fill=Species))+
      fill_colour_TL_4+
      theme
legend<-get_legend(p3)

p3<-ggplot(data=databis2[databis2$ea==as.character(ea10) & databis2$ma==as.character(ma10) & databis2$model=="4",])+
      geom_rect(aes(xmin=x-0.5,xmax=x+0.5,ymin=1e-3,ymax=Biomass,fill=Species))+
      facet_grid(ma~ea, labeller=label_parsed)+
      fill_colour_TL_4+
      y_axis_log10_short+
      coord_flip()+
      theme+theme(legend.position = "none")+
      xlab('Trophic level')+
      ylab("Biomass")

p4<-ggplot(data=databis)+
      geom_line(aes(S0_pred,Correlation,colour=Species,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_3+
      scale_linetype_manual(values=c("dotted","solid","twodash"),
                            guide = guide_legend(reverse = TRUE),
                            name='Chain\nlength')+
      theme+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_S0)+
      ylab(label_correlation)+
      ggtitle(expression("Dispersal of prey ("*italic(d["i"])*"=10"^"3"*")"))

### CV - DEMOGRAPHIC PERTURBATIONS ####
# parameters
nSpecies=4
nCommunity=1
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of species 1
                         disp=list(c(0,0,0,0)), # no dispersal
                         d=0,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=1,
                         C0=0,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="pert_1")
params_data<-rbind(params_data,expand.grid(pert=list(list(c(2,1))), # perturbation of species 2
                                           disp=list(c(0,0,0,0)), # no dispersal
                                           d=0,
                                           S0_self=1,
                                           S0_prey=1,
                                           S0_pred=1,
                                           C0=0,
                                           C0_self=0,
                                           C0_prey=0,
                                           C0_pred=0,
                                           correction=correction,
                                           model="pert_2"))
params_data<-rbind(params_data,expand.grid(pert=list(list(c(3,1))), # perturbation of species 3
                                           disp=list(c(0,0,0,0)), # no dispersal
                                           d=0,
                                           S0_self=1,
                                           S0_prey=1,
                                           S0_pred=1,
                                           C0=0,
                                           C0_self=0,
                                           C0_prey=0,
                                           C0_pred=0,
                                           correction=correction,
                                           model="pert_3"))
params_data<-rbind(params_data,expand.grid(pert=list(list(c(4,1))), # perturbation of species 4
                                           disp=list(c(0,0,0,0)), # no dispersal
                                           d=0,
                                           S0_self=1,
                                           S0_prey=1,
                                           S0_pred=1,
                                           C0=0,
                                           C0_self=0,
                                           C0_prey=0,
                                           C0_pred=0,
                                           correction=correction,
                                           model="pert_4"))
params_data<-merge(params_data_original,params_data)
params_data<-params_data[which(params_data$ea%in%as.character(c(ea01,ea10))),]
params_data<-params_data[which(params_data$ma%in%as.character(c(ma01,ma10))),]
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-data_V[,-which(names(data_V)%in%names(params_data)[-which(names(params_data)%in%c("ea","ma","model"))])]
databis<-table_for_matrix(databis,3)
databis<-databis[databis$species_1==databis$species_2,]
databis$species_2<-NULL
names(databis)[names(databis)=="species_1"]="species"
databis_2<-data_B[,-which(names(data_B)%in%names(params_data)[-which(names(params_data)%in%c("ea","ma","model"))])]
databis_2<-table_for_biomass(databis_2,3)
databis<-merge(databis,databis_2,by=c("ea","ma","model","species"))
rm(databis_2)
databis$CV<-sqrt(databis$value)/databis$biomass
levels(databis$model)<-c(1:nSpecies)
databis$model<-as.numeric(databis$model)

p5<-ggplot(data=databis)+
      geom_vline(aes(xintercept=model, colour=as.factor(model)),linetype='dashed',show.legend = FALSE)+
      geom_line(aes(model,CV,colour=species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      y_axis_log10_short+
      theme+
      #ggtitle("Single patch")+
      xlab("Perturbed trophic level")+
      ylab("Coefficient of variation (CV)")

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 1.6)) +
  draw_plot(p1, 0, 1, 0.6, 0.6)+
  draw_plot(p2, 0.6, 1, 0.6, 0.6)+
  draw_plot(p3, 1.2, 1, 0.6, 0.6)+
  draw_plot(legend, 1.85, 1.05, 0.05, 0.5)+
  draw_plot(p4, 0, 0, 1, 1)+
  draw_plot(p5, 1, 0, 1, 1)+
  draw_plot_label(c("A","B","C"), c(0,0,1), c(1.6,1,1), size = 30)
ggsave(paste(path_figure_results,"figure_CV_top_down.pdf",sep=""),graph, width = 18, height = 13, device = cairo_pdf)

########################## ----
# REALISTIC CASE - ENVIRONMENTAL PERTURBATION #### ----
### SIMULATIONS - d - CORRELATION - PASSIVE DISPERSAL - CORRELATED ENVIRONMENTAL PERTURBATIONS ####
# parameters
nSpecies=4
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(2,1),c(3,1),c(4,1))), # perturbation of predators in patch 1
                         disp=list(c(1,1,1,1)), # dispersal of prey
                         d=d,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=1,
                         C0=1,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="disp_passive")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
params_data$z=1 # environmental perturbation
VE<-matrix(sigma^2,nSpecies*nCommunity,nSpecies*nCommunity) # correlated perturbations
params_data$VE<-list(VE)
params_data<-params_data[params_data$ea==as.character(ea01) & params_data$ma==as.character(ma10),]

# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","C_11_12",'C_21_22','C_31_32','C_41_42'))],
              id.vars = c("ea","ma","d"),
              variable.name = "species",
              value.name = "correlation")

p1<-ggplot(data=databis)+
      geom_line(aes(d,correlation,colour=species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle("Passive dispersal")

# BIOMASS CV
databis<-data_V[,-which(names(data_V)%in%names(params_data)[-which(names(params_data)%in%c("ea","ma","d"))])]
databis<-table_for_matrix(databis,3)
databis<-databis[databis$species_1==databis$species_2 & databis$community_1==databis$community_2,]
databis$species_2<-NULL
databis$community_2<-NULL
names(databis)[names(databis)=="species_1"]="species"
names(databis)[names(databis)=="community_1"]="community"
databis_2<-data_B[,-which(names(data_B)%in%names(params_data)[-which(names(params_data)%in%c("ea","ma","d"))])]
databis_2<-table_for_biomass(databis_2,3)
databis<-merge(databis,databis_2,by=c("ea","ma","d","species","community"))
rm(databis_2)
databis$CV<-sqrt(databis$value)/databis$biomass
databis$species<-as.factor(databis$species)
databis$community<-as.factor(databis$community)

p3<-ggplot(data=databis)+
      geom_line(aes(d,CV,colour=species,linetype=community),size=1.1)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      patch_line+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      y_axis_log10_short+
      xlab(label_dispersal)+
      ylab(label_CV)

### SIMULATIONS - d - CORRELATION - EVERYTHING - CONTRIBUTION DEPENDING ON GROWTH RATE - CORRELATED ENVIRONMENTAL PERTURBATIONS ####
# parameters
nSpecies=4
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(2,1),c(3,1),c(4,1))), # perturbation of predators in patch 1
                         disp=list(c(1,1,1,1)), # dispersal of prey
                         d=d,
                         S0_self=1e3,
                         S0_prey=1e-3,
                         S0_pred=1e3,
                         C0=0,
                         C0_self=1,
                         C0_prey=1,
                         C0_pred=1,
                         correction=correction,
                         model="disp_prey")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
params_data$z=1 # environmental perturbation
VE<-matrix(sigma^2,nSpecies*nCommunity,nSpecies*nCommunity) # correlated perturbations
params_data$VE<-list(VE)
params_data$weight=1
params_data<-params_data[params_data$ea==as.character(ea01) & params_data$ma==as.character(ma10),]
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","C_11_12",'C_21_22','C_31_32','C_41_42'))],
              id.vars = c("ea","ma","d"),
              variable.name = "species",
              value.name = "correlation")

p2<-ggplot(data=databis)+
      geom_line(aes(d,correlation,colour=species),size=1.5)+
      corr_colour_TL_4+
      theme
legend_1<-get_legend(p2)

p2<-ggplot(data=databis)+
      geom_line(aes(d,correlation,colour=species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle("Weighted dependencies")

# BIOMASS CV
databis<-data_V[,-which(names(data_V)%in%names(params_data)[-which(names(params_data)%in%c("ea","ma","d"))])]
databis<-table_for_matrix(databis,3)
databis<-databis[databis$species_1==databis$species_2 & databis$community_1==databis$community_2,]
databis$species_2<-NULL
databis$community_2<-NULL
names(databis)[names(databis)=="species_1"]="species"
names(databis)[names(databis)=="community_1"]="community"
databis_2<-data_B[,-which(names(data_B)%in%names(params_data)[-which(names(params_data)%in%c("ea","ma","d"))])]
databis_2<-table_for_biomass(databis_2,3)
databis<-merge(databis,databis_2,by=c("ea","ma","d","species","community"))
rm(databis_2)
databis$CV<-sqrt(databis$value)/databis$biomass
databis$species<-as.factor(databis$species)
databis$community<-as.factor(databis$community)

p4<-ggplot(data=databis)+
      geom_line(aes(d,CV,colour=species,linetype=community),size=1.1)+
      corr_colour_TL_4+
      patch_line+
      theme
legend_2<-get_legend(p4)

p4<-ggplot(data=databis)+
      geom_line(aes(d,CV,colour=species,linetype=community),size=1.1)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      patch_line+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      y_axis_log10_short+
      xlab(label_dispersal)+
      ylab(label_CV)

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(p4, 1, 0, 1, 1)+
  draw_plot(legend_1, 2.05, 1.25, 0.05, 0.5)+
  draw_plot(legend_2, 2.05, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 30)
ggsave(paste(path_figure_results,"figure_env_everything.pdf",sep=""),graph, width = 14, height = 12, device = cairo_pdf)

########################## ----
# RELATIVE IMPORTANCE OF DISPERSAL ####
nSpecies=4
nCommunity=1
params_data<-expand.grid(d=d,
                         S0_self=1e3,
                         S0_prey=1e-3,
                         S0_pred=1e3,
                         C0=1,
                         C0_self=1,
                         C0_prey=1,
                         C0_pred=1,
                         correction=correction)
params_data<-merge(params_data_original,params_data)
params_data<-params_data[params_data$ea==as.character(ea01) & params_data$ma==as.character(ma10),]
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% get_dispersal_importance(params_data,nSpecies,i)
stopCluster(cl)
data<-params_data
data$'1'=0
data$'2'=0
data$'3'=0
data$'4'=0
for (i in 1:dim(params_data)[1]){
  data[i,which(names(data)%in%as.character(1:nSpecies))]<-unlist(results[[i]])
}

databis<-melt(data[,which(names(data)%in%c("ea","ma","d","1",'2',"3","4"))],
              id.vars = c("ea","ma","d"),
              variable.name = "species",
              value.name = "M")

p1<-ggplot(data=databis)+
      geom_line(aes(d,M,colour=species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      x_axis_log10_short+
      theme+
      xlab(label_dispersal)+
      ylab(label_balance)

# WEIGHT DEPENDENCY ####
# parameters
nSpecies=4
nCommunity=1
params_data<-expand.grid(pert=list(list(c(1,1),c(2,1),c(3,1),c(4,1))), # perturbation of species 1
                         disp=list(c(0,0,0,0)), # no dispersal
                         d=0,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=1,
                         C0=0,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="all")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
params_data$z=1 # environmental perturbation
VE<-matrix(sigma^2,nSpecies*nCommunity,nSpecies*nCommunity) # correlated perturbations
params_data$VE<-list(VE)
params_data<-params_data[params_data$ea==as.character(ea01) & params_data$ma==as.character(ma10),]
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-NULL
for(i in 1:dim(data_B)[1]){
  weight<-get_weight(as.numeric(data_B[i,which(names(data_B)%in%c("B_11","B_21","B_31","B_41"))]),
                     data_B$e[i],data_B$a[i],data_B$m[i],
                     1,1,1,
                     nSpecies)
  weight$e=data_B$e[i]
  weight$a=data_B$a[i]
  weight$m=data_B$m[i]
  databis<-rbind(databis,weight)
}
databis<-merge(databis,data_B[,which(names(data_B)%in%c("e","a","m","ea","ma"))],
               by=c("e","a","m"))

p2<-ggplot(data=databis)+
      geom_col(aes(species,weight,fill=dependency),position="fill")+
      facet_grid(ma~ea, labeller=label_parsed)+
      fill_colour_weight+
      theme+
      xlab('Trophic level')+
      ylab('Weight of dispersal components')

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure_results,"figure_env_M_weight.pdf",sep=""),graph, width = 16, height = 6, device = cairo_pdf)

########################## ----
# APPENDIX ############### ----
########################## ----
# PASSIVE DISPERSAL ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         d=d,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=1,
                         C0=1,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="disp_pred")
params_data<-rbind(params_data,expand.grid(pert=list(list(c(2,1))), # perturbation of predators in patch 1
                                           disp=list(c(1,0)), # dispersal of prey
                                           d=d,
                                           S0_self=1,
                                           S0_prey=1,
                                           S0_pred=1,
                                           C0=1,
                                           C0_self=0,
                                           C0_prey=0,
                                           C0_pred=0,
                                           correction=correction,
                                           model="disp_prey"))
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","model","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","model"),
              variable.name = "species",
              value.name = "correlation")

p1<-ggplot(data=databis[databis$model=="disp_pred",])+
      geom_line(aes(d,correlation,colour=species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle("Passive dispersal of predators")

p2<-ggplot(data=databis[databis$model=="disp_prey",])+
      geom_line(aes(d,correlation,colour=species),size=1.5)+
      corr_colour_TL_2+
      theme
legend<-get_legend(p2)

p2<-ggplot(data=databis[databis$model=="disp_prey",])+
      geom_line(aes(d,correlation,colour=species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle("Passive dispersal of prey")

p1_img<-image_read_pdf(paste(path_figure_supp,"schema_little_pert_TL1_disp_TL2.pdf",sep=""))
p2_img<-image_read_pdf(paste(path_figure_supp,"schema_little_pert_TL2_disp_TL1.pdf",sep=""))

graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 1.27)) +
  draw_image(p1_img, 0.25, 1.02, 0.4, 0.25)+
  draw_image(p2_img, 1.25, 1.02, 0.4, 0.25)+
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure_supp,"supp_simple.pdf",sep=""),graph, width = 17, height = 10, device = cairo_pdf)

### CORRELATION MATRIX ####
databis<-data_C[data_C$model=="disp_pred"
                & data_C$d==max(data_C$d)
                & data_C$ma==as.character(ma10)
                & data_C$ea==as.character(ea01), c(which(names(data_C)%in%c("ea","ma")),(which(names(data_C)%in%"VE")+1):dim(data_C)[2])]
databis<-melt(databis,
              id.vars = c("ea","ma"),
              variable.name = "species",
              value.name = "correlation")
databis<-databis %>% separate(species,c(NA,"species_1","species_2"),sep="_")
databis$species_1<-as.factor(databis$species_1)
databis$species_2<-as.factor(databis$species_2)
databis$species_1<-factor(databis$species_1,levels(databis$species_1)[c(1,3,2,4)])
databis$species_2<-factor(databis$species_2,levels(databis$species_2)[c(1,3,2,4)])
levels(databis$species_1)<-c("1","2"," 1"," 2")
levels(databis$species_2)<-c("1","2"," 1"," 2")

graph<-ggplot(data=databis)+
          geom_raster(aes(species_2,species_1,fill=correlation))+
          geom_vline(xintercept=2.5,linetype='dashed',size=1.5,colour='white')+
          geom_hline(yintercept=2.5,linetype='dashed',size=1.5,colour='white')+
          # annotate('text', colour='white',size=6, family="bold",
          #          x=c(1.5,3.5,1.5,3.5),
          #          y=c(1.5,3.5,3.5,1.5),
          #          label=c('intra\n#1','intra\n#2','inter','inter'))+
          #facet_grid(ma~ea, labeller=label_parsed)+
          corr_colour_grad+
          theme_matrix+theme(axis.ticks = element_blank(),
                             axis.text = element_blank(),
                             axis.line = element_blank())

ggsave(paste(path_figure_supp,"supp_correlation_matrix.pdf",sep=""),graph, width = 5, height = 3.5, device = cairo_pdf)

### TIME SERIES PULSE ####
#parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         d=1,
                         S0_self=1,
                         S0_prey=1e-3,
                         S0_pred=1,
                         C0=0,
                         C0_self=0,
                         C0_prey=1,
                         C0_pred=0,
                         correction=correction,
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
params_data<-params_data[params_data$ma==as.character(ma10) & params_data$ea==as.character(ea01),]
time<-seq(0,0.10,by=1e-3)
# simulations
data_TS<-ODE_solve(time, params_data, nSpecies, nCommunity, 1)
#data_TS<-merge(data_TS[seq(1,dim(data_TS)[1],by=1e1),],params_data)
data_TS<-merge(data_TS,params_data)

databis<-melt(data_TS[,-which(names(data_TS)%in%names(params_data[-which(names(params_data)%in%c("ma","ea","time"))]))],
              id.vars = c("ma","ea","time"),
              variable.name = "species",
              value.name = "biomass")
levels(databis$species)<-c("1_1","2_1","1_2","2_2")
databis<-databis %>% separate(species,c("species","community"),sep="_")

graph<-ggplot(data=databis)+
          geom_line(aes(time,biomass,colour=species,linetype=community),size=1.5)+
          facet_grid(ma~ea, labeller=label_parsed, scales = "free")+
          corr_colour_TL_2+
          patch_line+
          # corr_colour_TL_2_TS+
          # corr_line_TL_2_TS+
          theme+#theme(legend.position = "none")+
          y_axis_log10+
          xlab("Time")+
          ylab("Biomass")+
          ggtitle(expression(atop("Dispersal of predators","depending on prey")))

ggsave(paste(path_figure_supp,"supp_time_series.pdf",sep=""),graph, width = 7, height = 6, device = cairo_pdf)

########################## ----
# DISPERSAL DEPENDING ON SELF INTERACTIONS ####
### SIMULATIONS - d - CORRELATION ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         d=d,
                         S0_self=c(1e-3,1e3),
                         S0_prey=1,
                         S0_pred=1,
                         C0=0,
                         C0_self=1,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="disp_pred")
params_data<-rbind(params_data,expand.grid(pert=list(list(c(2,1))), # perturbation of predators in patch 1
                                           disp=list(c(1,0)), # dispersal of prey
                                           d=d,
                                           S0_self=c(1e-3,1e3),
                                           S0_prey=1,
                                           S0_pred=1,
                                           C0=0,
                                           C0_self=1,
                                           C0_prey=0,
                                           C0_pred=0,
                                           correction=correction,
                                           model="disp_prey"))
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_self","model","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","S0_self","model"),
              variable.name = "Species",
              value.name = "Correlation")

p1<-ggplot(data=databis[databis$S0_self==1e-3 & databis$model=="disp_pred",])+
      geom_line(aes(d,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle(expression("Dispersal of predators ("*italic(S["0"])*"=10"^"-3"*")"))

# p2<-ggplot(data=databis[databis$S0_self==1e-3 & databis$model=="disp_prey",])+
#       geom_line(aes(d,Correlation,colour=Species),size=1.5)+
#       facet_grid(ma~ea, labeller=label_parsed)+
#       corr_colour_TL_2+
#       theme+theme(legend.position = "none")+
#       x_axis_log10_short+
#       ylim(-1,1)+
#       xlab(label_dispersal)+
#       ylab(label_correlation)+
#       ggtitle("Dispersal of prey - S0=1e-3")

# p3<-ggplot(data=databis[databis$S0_self==1e3 & databis$model=="disp_pred",])+
#       geom_line(aes(d,Correlation,colour=Species),size=1.5)+
#       facet_grid(ma~ea, labeller=label_parsed)+
#       corr_colour_TL_2+
#       theme+theme(legend.position = "none")+
#       x_axis_log10_short+
#       ylim(-1,1)+
#       xlab(label_dispersal)+
#       ylab(label_correlation)+
#       ggtitle("Dispersal of predators - S0=1e3")

p2<-ggplot(data=databis[databis$S0_self==1e3 & databis$model=="disp_prey",])+
      geom_line(aes(d,Correlation,colour=Species),size=1.5)+
      corr_colour_TL_2+
      theme
legend<-get_legend(p2)

p2<-ggplot(data=databis[databis$S0_self==1e3 & databis$model=="disp_prey",])+
      geom_line(aes(d,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle(expression(paste("Dispersal of prey (",italic(S["0"]),"=10"^"3",")")))

# graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 2)) +
#   draw_plot(p1, 0, 1, 1, 1)+
#   draw_plot(p2, 1, 1, 1, 1)+
#   draw_plot(p3, 0, 0, 1, 1)+
#   draw_plot(p4, 1, 0, 1, 1)+
#   draw_plot(legend, 2.05, 0.25, 0.05, 0.5)+
#   draw_plot(legend, 2.05, 1.25, 0.05, 0.5)+
#   draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 30)
# ggsave(paste(path_figure,"figure_self_d.pdf",sep=""),graph, width = 17, height = 16, device = cairo_pdf)

### SIMULATIONS - S0 - CORRELATION ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         d=1e3,
                         S0_self=S0,
                         S0_prey=1,
                         S0_pred=1,
                         C0=0,
                         C0_self=1,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="disp_pred")
params_data<-rbind(params_data,expand.grid(pert=list(list(c(2,1))), # perturbation of predators in patch 1
                                           disp=list(c(1,0)), # dispersal of prey
                                           d=1e3,
                                           S0_self=S0,
                                           S0_prey=1,
                                           S0_pred=1,
                                           C0=0,
                                           C0_self=1,
                                           C0_prey=0,
                                           C0_pred=0,
                                           correction=correction,
                                           model="disp_prey"))
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_self","model","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","S0_self","model"),
              variable.name = "Species",
              value.name = "Correlation")

p3<-ggplot(data=databis[databis$d==1e3 & databis$model=="disp_pred",])+
      geom_line(aes(S0_self,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_S0)+
      ylab(label_correlation)+
      ggtitle(expression("Dispersal of predators ("*italic(d["i"])*"=10"^"3"*")"))

# p2<-ggplot(data=databis[databis$d==1e3 & databis$model=="disp_prey",])+
#       geom_line(aes(S0_self,Correlation,colour=Species),size=1.5)+
#       corr_colour_TL_2
# legend<-get_legend(p2)

p4<-ggplot(data=databis[databis$d==1e3 & databis$model=="disp_prey",])+
      geom_line(aes(S0_self,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_S0)+
      ylab(label_correlation)+
      ggtitle(expression("Dispersal of prey ("*italic(d["i"])*"=10"^"3"*")"))

# graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 1)) +
#   draw_plot(p1, 0, 0, 1, 1)+
#   draw_plot(p2, 1, 0, 1, 1)+
#   draw_plot(legend, 2.05, 0.25, 0.05, 0.5)+
#   draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
# ggsave(paste(path_figure,"figure_self_S0.pdf",sep=""),graph, width = 17, height = 8, device = cairo_pdf)

### FINAL GRAPH ####

graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(p4, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 0.25, 0.05, 0.5)+
  draw_plot(legend, 2.05, 1.25, 0.05, 0.5)+
  draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 30)
ggsave(paste(path_figure_supp,"supp_self.pdf",sep=""),graph, width = 17, height = 16, device = cairo_pdf)

########################## ----
# DENSITY DEPENDENT DISPERSAL FULL PARAMETERS #### ----
# DISPERSAL OF PRED DEPENDING ON PREY #### ----
### SIMULATIONS - d - CORRELATION ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         d=d,
                         S0_self=1,
                         S0_prey=c(1e-3,1e3),
                         S0_pred=1,
                         C0=0,
                         C0_self=0,
                         C0_prey=1,
                         C0_pred=0,
                         correction=correction,
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_prey","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","S0_prey"),
              variable.name = "Species",
              value.name = "Correlation")

p1<-ggplot(data=databis[databis$S0_prey==1e-3,])+
      geom_vline(xintercept=1e3,linetype='dashed')+
      geom_line(aes(d,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle(expression(atop("Dispersal of predators","depending on prey ("*italic(S["0"])*"=10"^"-3"*")")))


p3<-ggplot(data=databis[databis$S0_prey==1e3,])+
      geom_vline(xintercept=1e3,linetype='dashed')+
      geom_line(aes(d,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle(expression(atop("Dispersal of predators","depending on prey ("*italic(S["0"])*"=10"^"3"*")")))

### SIMULATIONS - S0 - CORRELATION ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         d=1e3,
                         S0_self=1,
                         S0_prey=S0,
                         S0_pred=1,
                         C0=0,
                         C0_self=0,
                         C0_prey=1,
                         C0_pred=0,
                         correction=correction,
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_prey","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","S0_prey"),
              variable.name = "Species",
              value.name = "Correlation")

p5<-ggplot(data=databis[databis$d==1e3,])+
      geom_line(aes(S0_prey,Correlation,colour=Species),size=1.5)+
      corr_colour_TL_2+
      theme
legend<-get_legend(p5)

p5<-ggplot(data=databis[databis$d==1e3,])+
      geom_vline(xintercept=1e-3,linetype='dashed')+
      geom_vline(xintercept=1e3,linetype='dashed')+
      geom_line(aes(S0_prey,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_S0)+
      ylab(label_correlation)+
      ggtitle(expression(atop("Dispersal of predators","depending on prey ("*italic(d["i"])*"=10"^"3"*")")))

# DISPERSAL OF PREY DEPENDING ON PRED #### ----
### SIMULATIONS - d - CORRELATION ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(2,1))), # perturbation of predators in patch 1
                         disp=list(c(1,0)), # dispersal of prey
                         d=d,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=c(1e-3,1e3),
                         C0=0,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=1,
                         correction=correction,
                         model="disp_prey")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_pred","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","S0_pred"),
              variable.name = "Species",
              value.name = "Correlation")

p2<-ggplot(data=databis[databis$S0_pred==1e-3,])+
      geom_vline(xintercept=1e3,linetype='dashed')+
      geom_line(aes(d,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle(expression(atop("Dispersal of prey","depending on predators ("*italic(S["0"])*"=10"^"-3"*")")))

p4<-ggplot(data=databis[databis$S0_pred==1e3,])+
      geom_vline(xintercept=1e3,linetype='dashed')+
      geom_line(aes(d,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle(expression(atop("Dispersal of prey","depending on predators ("*italic(S["0"])*"=10"^"3"*")")))

### SIMULATIONS - S0 - CORRELATION ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(2,1))), # perturbation of predators in patch 1
                         disp=list(c(1,0)), # dispersal of prey
                         d=1e3,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=S0,
                         C0=0,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=1,
                         correction=correction,
                         model="disp_prey")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_pred","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","S0_pred"),
              variable.name = "Species",
              value.name = "Correlation")

p6<-ggplot(data=databis[databis$d==1e3,])+
      geom_vline(xintercept=1e-3,linetype='dashed')+
      geom_vline(xintercept=1e3,linetype='dashed')+
      geom_line(aes(S0_pred,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_S0)+
      ylab(label_correlation)+
      ggtitle(expression(atop("Dispersal of prey","depending on predators ("*italic(d["i"])*"=10"^"3"*")")))

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 3)) +
  draw_plot(p1, 0, 2, 1, 1)+
  draw_plot(p2, 1, 2, 1, 1)+
  draw_plot(p3, 0, 1, 1, 1)+
  draw_plot(p4, 1, 1, 1, 1)+
  draw_plot(p5, 0, 0, 1, 1)+
  draw_plot(p6, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 2.25, 0.05, 0.5)+
  draw_plot(legend, 2.05, 1.25, 0.05, 0.5)+
  draw_plot(legend, 2.05, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B","C","D","E","F"), c(0,1,0,1,0,1), c(3,3,2,2,1,1), size = 30)
ggsave(paste(path_figure_supp,"supp_prey_pred.pdf",sep=""),graph, width = 17, height = 24, device = cairo_pdf)

########################## ----
# S0 AND REGIONAL CV #### ----
# DISPERSAL OF PRED DEPENDING ON PREY ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         d=1e3,
                         S0_self=1,
                         S0_prey=S0,
                         S0_pred=1,
                         C0=0,
                         C0_self=0,
                         C0_prey=1,
                         C0_pred=0,
                         correction=correction,
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-table_for_matrix(data_V[,c(which(names(data_V)%in%c("ea","ma","S0_prey")),(which(names(data_V)%in%"VE")+1):dim(data_V)[2])],3)
databis<-databis[databis$species_1==databis$species_2,]
databis$species_2<-NULL
databis<-aggregate(databis$value, by=list(ea=databis$ea,
                                          ma=databis$ma,
                                          S0_prey=databis$S0_prey,
                                          species=databis$species_1), FUN=sum)
names(databis)[names(databis)=="x"]="variance"

databis_2<-table_for_biomass(data_B[,c(which(names(data_B)%in%c("ea","ma","S0_prey")),(which(names(data_B)%in%"VE")+1):dim(data_B)[2])],3)
databis_2<-aggregate(databis_2$biomass, by=list(ea=databis_2$ea,
                                                ma=databis_2$ma,
                                                S0_prey=databis_2$S0_prey,
                                                species=databis_2$species), FUN=sum)
names(databis_2)[names(databis_2)=="x"]="biomass"

databis<-merge(databis,databis_2,by=c("ea","ma","S0_prey","species"))
rm(databis_2)
databis$CV<-sqrt(databis$variance)/databis$biomass

p1<-ggplot(data=databis)+
      geom_line(aes(S0_prey,CV,colour=species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      xlab(label_S0)+
      ylab(label_CV)+
      ggtitle(expression(atop("Dispersal of predators","depending on prey ("*italic(d["i"])*"=10"^"3"*")")))

# DISPERSAL OF PREY DEPENDING ON PRED ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(2,1))), # perturbation of predators in patch 1
                         disp=list(c(1,0)), # dispersal of prey
                         d=1e3,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=S0,
                         C0=0,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=1,
                         correction=correction,
                         model="disp_prey")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-table_for_matrix(data_V[,c(which(names(data_V)%in%c("ea","ma","S0_pred")),(which(names(data_V)%in%"VE")+1):dim(data_V)[2])],3)
databis<-databis[databis$species_1==databis$species_2,]
databis$species_2<-NULL
databis<-aggregate(databis$value, by=list(ea=databis$ea,
                                          ma=databis$ma,
                                          S0_pred=databis$S0_pred,
                                          species=databis$species_1), FUN=sum)
names(databis)[names(databis)=="x"]="variance"

databis_2<-table_for_biomass(data_B[,c(which(names(data_B)%in%c("ea","ma","S0_pred")),(which(names(data_B)%in%"VE")+1):dim(data_B)[2])],3)
databis_2<-aggregate(databis_2$biomass, by=list(ea=databis_2$ea,
                                                ma=databis_2$ma,
                                                S0_pred=databis_2$S0_pred,
                                                species=databis_2$species), FUN=sum)
names(databis_2)[names(databis_2)=="x"]="biomass"

databis<-merge(databis,databis_2,by=c("ea","ma","S0_pred","species"))
rm(databis_2)
databis$CV<-sqrt(databis$variance)/databis$biomass

p2<-ggplot(data=databis)+
  geom_line(aes(S0_pred,CV,colour=species),size=1.5)+
  corr_colour_TL_4+
  theme
legend<-get_legend(p2)

p2<-ggplot(data=databis)+
      geom_line(aes(S0_pred,CV,colour=species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      xlab(label_S0)+
      ylab(label_CV)+
      ggtitle(expression(atop("Dispersal of prey depending","on predators ("*italic(d["i"])*"=10"^"3"*")")))

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure_supp,"supp_prey_pred_CV.pdf",sep=""),graph, width = 17, height = 8, device = cairo_pdf)

########################## ----
# TOP-DOWN CONTROL SYMMETRIC #### ----
### SIMULATIONS - d - CORRELATION ####
# parameters
nSpecies=3
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1,0)), # dispersal of predators
                         d=d,
                         S0_self=1,
                         S0_prey=c(1e-3,1e3),
                         S0_pred=1,
                         C0=0,
                         C0_self=0,
                         C0_prey=1,
                         C0_pred=0,
                         correction=correction,
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_prey","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","S0_prey"),
              variable.name = "Species",
              value.name = "Correlation")

p1<-ggplot(data=databis[databis$S0_prey==1e-3,])+
      geom_vline(xintercept=1e3,linetype='dashed')+
      geom_line(aes(d,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle(expression(atop("Dispersal of predators","depending on prey ("*italic(S["0"])*"=10"^"-3"*")")))

### SIMULATIONS - S0 - CORRELATION - TOP-DOWN CONTROL - CHAIN LENGTH ####
# parameters - 2 species
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         d=1e3,
                         S0_self=1,
                         S0_prey=S0,
                         S0_pred=1,
                         C0=0,
                         C0_self=0,
                         C0_prey=1,
                         C0_pred=0,
                         correction=correction,
                         model="2")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations - 2 species
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_prey","model","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","S0_prey","model"),
              variable.name = "Species",
              value.name = "Correlation")

data_B$B_31=NA
data_B$B_41=NA
databis2<-melt(data_B[,which(names(data_B)%in%c("ea","ma","d","S0_prey","model","B_11",'B_21',"B_31","B_41"))],
               id.vars = c("ea","ma","d","S0_prey","model"),
               variable.name = "Species",
               value.name = "Biomass")

# parameters - 3 species
nSpecies=3
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1,0)), # dispersal of predators
                         d=1e3,
                         S0_self=1,
                         S0_prey=S0,
                         S0_pred=1,
                         C0=0,
                         C0_self=0,
                         C0_prey=1,
                         C0_pred=0,
                         correction=correction,
                         model="3")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations - 3 species
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-rbind(databis,melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_prey","model","C_11_12",'C_21_22'))],
                            id.vars = c("ea","ma","d","S0_prey","model"),
                            variable.name = "Species",
                            value.name = "Correlation"))

data_B$B_41=NA
databis2<-rbind(databis2,melt(data_B[,which(names(data_B)%in%c("ea","ma","d","S0_prey","model","B_11",'B_21',"B_31","B_41"))],
                              id.vars = c("ea","ma","d","S0_prey","model"),
                              variable.name = "Species",
                              value.name = "Biomass"))

# parameters - 4 species
nSpecies=4
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1,0,0)), # dispersal of predators
                         d=1e3,
                         S0_self=1,
                         S0_prey=S0,
                         S0_pred=1,
                         C0=0,
                         C0_self=0,
                         C0_prey=1,
                         C0_pred=0,
                         correction=correction,
                         model="4")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations - 4 species
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-rbind(databis,databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","S0_prey","model","C_11_12",'C_21_22'))],
                                     id.vars = c("ea","ma","d","S0_prey","model"),
                                     variable.name = "Species",
                                     value.name = "Correlation"))

databis2<-rbind(databis2,melt(data_B[,which(names(data_B)%in%c("ea","ma","d","S0_prey","model","B_11",'B_21',"B_31","B_41"))],
                              id.vars = c("ea","ma","d","S0_prey","model"),
                              variable.name = "Species",
                              value.name = "Biomass"))

p2<-ggplot(data=databis)+
      geom_line(aes(S0_prey,Correlation,colour=Species,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_3+
      scale_linetype_manual(values=c("dotted","solid","twodash"),
                            guide = guide_legend(reverse = TRUE),
                            name='Chain\nlength')+
      theme+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_S0)+
      ylab(label_correlation)+
      ggtitle(expression(atop("Dispersal of predators","depending on prey ("*italic(d["i"])*"=10"^"3"*")")))

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure_supp,"supp_CV_top_down.pdf",sep=""),graph, width = 16, height = 8, device = cairo_pdf)

########################## ----
# MULTIPLE TYPES OF DENSITY DEPENDENCY - 3 SPECIES #### ----
### SIMULATIONS - d - CORRELATION - PASSIVE DISPERSAL ####
# parameters
nSpecies=3
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(3,1))), # perturbation of PP and carnivores in patch 1
                         disp=list(c(0,1,0)), # dispersal of the herbivore
                         d=d,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=1,
                         C0=1,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="disp_herbi")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","C_11_12",'C_21_22','C_31_32'))],
              id.vars = c("ea","ma","d"),
              variable.name = "Species",
              value.name = "Correlation")

p1<-ggplot(data=databis)+
      geom_line(aes(d,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_3+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle("Passive dispersal of herbivores")

### SIMULATIONS - d - CORRELATION - DISPERSAL OF HERBIVORES DEPENDING ON SELF DEPENDENCY ####
# parameters
nSpecies=3
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(3,1))), # perturbation of PP and carnivores in patch 1
                         disp=list(c(0,1,0)), # dispersal of the herbivore
                         d=d,
                         S0_self=1e3,
                         S0_prey=1,
                         S0_pred=1,
                         C0=0,
                         C0_self=1,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="disp_herbi")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","C_11_12",'C_21_22','C_31_32'))],
              id.vars = c("ea","ma","d"),
              variable.name = "Species",
              value.name = "Correlation")

p2<-ggplot(data=databis)+
      geom_line(aes(d,Correlation,colour=Species),size=1.5)+
      corr_colour_TL_3
legend<-get_legend(p2)

p2<-ggplot(data=databis)+
      geom_line(aes(d,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_3+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle(expression("Dispersal of herbivores ("*italic(S["0"])*"=10"^"3"*")"))

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure_supp,"supp_simple_self.pdf",sep=""),graph, width = 17, height = 8, device = cairo_pdf)

########################## ----
# MULTIPLE TYPES OF DENSITY DEPENDENCY - 3 SPECIES #### ----
### REFERENCE VALUES BASED ON PASSIVE DISPERSAL ####
# parameters
nSpecies=3
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(3,1))), # perturbation of PP and carnivores in patch 1
                         disp=list(c(0,1,0)), # dispersal of the herbivore
                         d=d,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=1,
                         C0=1,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="passive")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

ref<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","model","C_11_12",'C_21_22','C_31_32'))],
          id.vars = c("ea","ma","d","model"),
          variable.name = "Species",
          value.name = "Correlation")

### SIMULATIONS - d - CORRELATION - DISPERSAL OF HERBIVORES DEPENDING ON PREY PLUS SELF DEPENDENCY - NOT WEIGHTED ####
# parameters
nSpecies=3
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(3,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1,0)), # dispersal of herbivoes
                         d=d,
                         S0_self=1e3,
                         S0_prey=1e-3,
                         S0_pred=1,
                         C0=0,
                         C0_self=1,
                         C0_prey=1,
                         C0_pred=0,
                         correction=correction,
                         model="density dep.")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","model","C_11_12",'C_21_22','C_31_32'))],
              id.vars = c("ea","ma","d","model"),
              variable.name = "Species",
              value.name = "Correlation")
databis<-rbind(databis,ref)

p1<-ggplot(data=databis)+
      geom_line(aes(d,Correlation,colour=Species,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_3+
      model_line+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle("Prey density dependent - not weighted")

### SIMULATIONS - d - CORRELATION - DISPERSAL OF HERBIVORES DEPENDING ON PREY PLUS SELF DEPENDENCY - WEIGHTED ####
# parameters
nSpecies=3
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(3,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1,0)), # dispersal of predators
                         d=d,
                         S0_self=1e3,
                         S0_prey=1e-3,
                         S0_pred=1,
                         C0=0,
                         C0_self=1,
                         C0_prey=1,
                         C0_pred=0,
                         correction=correction,
                         model="density dep.")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
params_data$weight=1
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","model","C_11_12",'C_21_22','C_31_32'))],
              id.vars = c("ea","ma","d","model"),
              variable.name = "Species",
              value.name = "Correlation")
databis<-rbind(databis,ref)

p2<-ggplot(data=databis)+
      geom_line(aes(d,Correlation,colour=Species,linetype=model),size=1.5)+
      corr_colour_TL_3+
      model_line+
      theme
legend<-get_legend(p2)

p2<-ggplot(data=databis)+
      geom_line(aes(d,Correlation,colour=Species,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_3+
      model_line+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle("Prey density dependent - weighted")

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 2.25), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 0.25, 0.15, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure_supp,"supp_prey_self.pdf",sep=""),graph, width = 17, height = 8, device = cairo_pdf)

########################## ----
# MULTIPLE TYPES OF DENSITY DEPENDENCY - 3 SPECIES #### ----
### REFERENCE VALUES BASED ON PASSIVE DISPERSAL ####
# parameters
nSpecies=3
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(3,1))), # perturbation of PP and carnivores in patch 1
                         disp=list(c(0,1,0)), # dispersal of the herbivore
                         d=d,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=1,
                         C0=1,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="passive")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

ref<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","model","C_11_12",'C_21_22','C_31_32'))],
          id.vars = c("ea","ma","d","model"),
          variable.name = "Species",
          value.name = "Correlation")

### SIMULATIONS - d - CORRELATION - DISPERSAL OF HERBIVORES DEPENDING ON PRED PLUS SELF DEPENDENCY - NOT WEIGHTED ####
# parameters
nSpecies=3
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(3,1))), # perturbation of predators in patch 1
                         disp=list(c(0,1,0)), # dispersal of prey
                         d=d,
                         S0_self=1e3,
                         S0_prey=1,
                         S0_pred=1e3,
                         C0=0,
                         C0_self=1,
                         C0_prey=0,
                         C0_pred=1,
                         correction=correction,
                         model="density dep.")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","model","C_11_12",'C_21_22','C_31_32'))],
              id.vars = c("ea","ma","d","model"),
              variable.name = "Species",
              value.name = "Correlation")
databis<-rbind(databis,ref)

p1<-ggplot(data=databis)+
      geom_line(aes(d,Correlation,colour=Species,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_3+
      model_line+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle("Predator density dependent - not weighted")

### SIMULATIONS - d - CORRELATION - DISPERSAL OF HERBIVORES DEPENDING ON PRED PLUS SELF DEPENDENCY - WEIGHTED ####
# parameters
nSpecies=3
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(3,1))), # perturbation of predators in patch 1
                         disp=list(c(0,1,0)), # dispersal of prey
                         d=d,
                         S0_self=1e3,
                         S0_prey=1,
                         S0_pred=1e3,
                         C0=0,
                         C0_self=1,
                         C0_prey=0,
                         C0_pred=1,
                         correction=correction,
                         model="density dep.")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
params_data$weight=1
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","model","C_11_12",'C_21_22','C_31_32'))],
              id.vars = c("ea","ma","d","model"),
              variable.name = "Species",
              value.name = "Correlation")
databis<-rbind(databis,ref)

p2<-ggplot(data=databis)+
      geom_line(aes(d,Correlation,colour=Species,linetype=model),size=1.5)+
      corr_colour_TL_3+
      model_line+
      theme
legend<-get_legend(p2)

p2<-ggplot(data=databis)+
      geom_line(aes(d,Correlation,colour=Species,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_3+
      model_line+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle("Predator density dependent - weighted")

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 2.25), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 0.25, 0.15, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure_supp,"supp_pred_self.pdf",sep=""),graph, width = 17, height = 8, device = cairo_pdf)

########################## ----
# MULTIPLE TYPES OF DENSITY DEPENDENCY - 3 SPECIES #### ----
### WEIGHT DEPENDENCY ####
# parameters
nSpecies=3
nCommunity=1
params_data<-expand.grid(pert=list(list(c(1,1),c(3,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1,0)), # dispersal of predators
                         d=0,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=1,
                         C0=0,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="prey dep.")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-NULL
for(i in 1:dim(data_B)[1]){
  weight<-get_weight(as.numeric(data_B[i,which(names(data_B)%in%c("B_11","B_21","B_31","B_41"))]),
                     data_B$e[i],data_B$a[i],data_B$m[i],
                     1,1,0,
                     nSpecies)
  weight$e=data_B$e[i]
  weight$a=data_B$a[i]
  weight$m=data_B$m[i]
  databis<-rbind(databis,weight)
}
databis<-merge(databis,data_B[,which(names(data_B)%in%c("e","a","m","ea","ma"))],
               by=c("e","a","m"))

p1<-ggplot(data=databis)+
  geom_col(aes(species,weight,fill=dependency),position="fill")+
  facet_grid(ma~ea, labeller=label_parsed)+
  fill_colour_weight+
  theme+theme(legend.position = "none")+
  xlab('Trophic level')+
  ylab('Weight of dispersal components')+
  ggtitle("Prey density dependent")

databis<-NULL
for(i in 1:dim(data_B)[1]){
  weight<-get_weight(as.numeric(data_B[i,which(names(data_B)%in%c("B_11","B_21","B_31","B_41"))]),
                     data_B$e[i],data_B$a[i],data_B$m[i],
                     1,0,1,
                     nSpecies)
  weight$e=data_B$e[i]
  weight$a=data_B$a[i]
  weight$m=data_B$m[i]
  databis<-rbind(databis,weight)
}
databis<-merge(databis,data_B[,which(names(data_B)%in%c("e","a","m","ea","ma"))],
               by=c("e","a","m"))

p2<-ggplot(data=databis)+
      geom_col(aes(species,weight,fill=dependency),position="fill")+
      fill_colour_weight+
      theme
legend<-get_legend(p2)

p2<-ggplot(data=databis)+
      geom_col(aes(species,weight,fill=dependency),position="fill")+
      facet_grid(ma~ea, labeller=label_parsed)+
      fill_colour_weight+
      theme+theme(legend.position = "none")+
      xlab('Trophic level')+
      ylab('Weight of dispersal components')+
      ggtitle("Predator density dependent")

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 2.25), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 0.25, 0.15, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure_supp,"supp_weight.pdf",sep=""),graph, width = 17, height = 8, device = cairo_pdf)

########################## ----
# MULTIPLE TYPES OF DENSITY DEPENDENCY - 3 SPECIES #### ----
### REFERENCE VALUES BASED ON PASSIVE DISPERSAL ####
# parameters
nSpecies=3
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(3,1))), # perturbation of PP and carnivores in patch 1
                         disp=list(c(0,1,0)), # dispersal of the herbivore
                         d=d,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=1,
                         C0=1,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="passive")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

ref<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","model","C_11_12",'C_21_22','C_31_32'))],
          id.vars = c("ea","ma","d","model"),
          variable.name = "Species",
          value.name = "Correlation")

### SIMULATIONS - d - CORRELATION - EVERYTHING - EQUAL CONTRIBUTION ####
# parameters
nSpecies=3
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(3,1))), # perturbation of predators in patch 1
                         disp=list(c(0,1,0)), # dispersal of prey
                         d=d,
                         S0_self=1e3,
                         S0_prey=1e-3,
                         S0_pred=1e3,
                         C0=0,
                         C0_self=1,
                         C0_prey=1,
                         C0_pred=1,
                         correction=correction,
                         model="density dep.")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","model","C_11_12",'C_21_22','C_31_32'))],
              id.vars = c("ea","ma","d","model"),
              variable.name = "Species",
              value.name = "Correlation")
databis<-rbind(databis,ref)

p1<-ggplot(data=databis)+
      geom_line(aes(d,Correlation,colour=Species,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_3+
      model_line+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle("Equally weighted dependencies")

### SIMULATIONS - d - CORRELATION - EVERYTHING - CONTRIBUTION DEPENDING ON GROWTH RATE ####
# parameters
nSpecies=3
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(3,1))), # perturbation of predators in patch 1
                         disp=list(c(0,1,0)), # dispersal of prey
                         d=d,
                         S0_self=1e3,
                         S0_prey=1e-3,
                         S0_pred=1e3,
                         C0=0,
                         C0_self=1,
                         C0_prey=1,
                         C0_pred=1,
                         correction=correction,
                         model="density dep.")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
params_data$weight=1
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","model","C_11_12",'C_21_22','C_31_32'))],
              id.vars = c("ea","ma","d","model"),
              variable.name = "Species",
              value.name = "Correlation")
databis<-rbind(databis,ref)

p2<-ggplot(data=databis)+
      geom_line(aes(d,Correlation,colour=Species,linetype=model),size=1.5)+
      corr_colour_TL_3+
      model_line+
      theme
legend<-get_legend(p2)

p2<-ggplot(data=databis)+
      geom_line(aes(d,Correlation,colour=Species,linetype=model),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_3+
      model_line+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle("Weighted dependencies")

### CV - DEMOGRAPHIC PERTURBATIONS - WEIGHT DEPENDENCY ####
# parameters
nSpecies=3
nCommunity=1
params_data<-expand.grid(pert=list(list(c(1,1),c(3,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1,0)), # dispersal of predators
                         d=0,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=1,
                         C0=0,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="density dep.")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-data_V[,-which(names(data_V)%in%names(params_data)[-which(names(params_data)%in%c("ea","ma","model"))])]
databis<-table_for_matrix(databis,3)
databis<-databis[databis$species_1==databis$species_2,]
databis$species_2<-NULL
names(databis)[names(databis)=="species_1"]="species"
databis_2<-data_B[,-which(names(data_B)%in%names(params_data)[-which(names(params_data)%in%c("ea","ma","model"))])]
databis_2<-table_for_biomass(databis_2,3)
databis<-merge(databis,databis_2,by=c("ea","ma","model","species"))
rm(databis_2)
databis$CV<-sqrt(databis$value)/databis$biomass
levels(databis$model)<-c(1:nSpecies)
databis$model<-as.numeric(databis$model)
databis$x<-as.numeric(databis$species)

p3<-ggplot(data=databis)+
      geom_rect(aes(xmin=x-0.5,xmax=x+0.5,ymin=1e-4,ymax=CV,fill=species))+
      facet_grid(ma~ea, labeller=label_parsed)+
      fill_colour_TL_3+
      y_axis_log10_short+
      theme+
      xlab('Trophic level')+
      ylab(label_CV)

databis<-NULL
for(i in 1:dim(data_B)[1]){
  weight<-get_weight(as.numeric(data_B[i,which(names(data_B)%in%c("B_11","B_21","B_31","B_41"))]),
                     data_B$e[i],data_B$a[i],data_B$m[i],
                     1,1,1,
                     nSpecies)
  weight$e=data_B$e[i]
  weight$a=data_B$a[i]
  weight$m=data_B$m[i]
  databis<-rbind(databis,weight)
}
databis<-merge(databis,data_B[,which(names(data_B)%in%c("e","a","m","ea","ma"))],
               by=c("e","a","m"))

p4<-ggplot(data=databis)+
      geom_col(aes(species,weight,fill=dependency),position="fill")+
      facet_grid(ma~ea, labeller=label_parsed)+
      fill_colour_weight+
      theme+
      xlab('Trophic level')+
      ylab('Weight of dispersal components')

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 2.25), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1.125, 1)+
  draw_plot(p4, 1.125, 0, 1.125, 1)+
  draw_plot(legend, 2.05, 1.25, 0.15, 0.5)+
  draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 30)
ggsave(paste(path_figure_supp,"supp_everything.pdf",sep=""),graph, width = 17, height = 16, device = cairo_pdf)

########################## ----
# REALISTIC CASE - ENVIRONMENTAL PERTURBATION #### ----
### SIMULATIONS - d - CORRELATION - PASSIVE DISPERSAL - CORRELATED ENVIRONMENTAL PERTURBATIONS ####
# parameters
nSpecies=4
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(2,1),c(3,1),c(4,1))), # perturbation of predators in patch 1
                         disp=list(c(1,1,1,1)), # dispersal of prey
                         d=d,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=1,
                         C0=1,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="disp_simple")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
params_data$z=1 # environmental perturbation
VE<-matrix(sigma^2,nSpecies*nCommunity,nSpecies*nCommunity) # correlated perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","C_11_12",'C_21_22','C_31_32','C_41_42'))],
              id.vars = c("ea","ma","d"),
              variable.name = "Species",
              value.name = "Correlation")

p1<-ggplot(data=databis)+
      geom_line(aes(d,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      theme+#theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle("Passive dispersal")

# BIOMASS CV
databis<-data_V[,-which(names(data_V)%in%names(params_data)[-which(names(params_data)%in%c("ea","ma","d"))])]
databis<-table_for_matrix(databis,3)
databis<-databis[databis$species_1==databis$species_2 & databis$community_1==databis$community_2,]
databis$species_2<-NULL
databis$community_2<-NULL
names(databis)[names(databis)=="species_1"]="species"
names(databis)[names(databis)=="community_1"]="community"
databis_2<-data_B[,-which(names(data_B)%in%names(params_data)[-which(names(params_data)%in%c("ea","ma","d"))])]
databis_2<-table_for_biomass(databis_2,3)
databis<-merge(databis,databis_2,by=c("ea","ma","d","species","community"))
rm(databis_2)
databis$CV<-sqrt(databis$value)/databis$biomass
databis$species<-as.factor(databis$species)
databis$community<-as.factor(databis$community)

p3<-ggplot(data=databis)+
      geom_line(aes(d,CV,colour=species,linetype=community),size=1.1)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      patch_line+
      theme+
      x_axis_log10_short+
      y_axis_log10_short+
      xlab(label_dispersal)+
      ylab(label_CV)

### SIMULATIONS - d - CORRELATION - EVERYTHING - CONTRIBUTION DEPENDING ON GROWTH RATE - CORRELATED ENVIRONMENTAL PERTURBATIONS ####
# parameters
nSpecies=4
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(2,1),c(3,1),c(4,1))), # perturbation of predators in patch 1
                         disp=list(c(1,1,1,1)), # dispersal of prey
                         d=d,
                         S0_self=1e3,
                         S0_prey=1e-3,
                         S0_pred=1e3,
                         C0=0,
                         C0_self=1,
                         C0_prey=1,
                         C0_pred=1,
                         correction=correction,
                         model="disp_prey")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
params_data$z=1 # environmental perturbation
VE<-matrix(sigma^2,nSpecies*nCommunity,nSpecies*nCommunity) # correlated perturbations
params_data$VE<-list(VE)
params_data$weight=1
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","C_11_12",'C_21_22','C_31_32','C_41_42'))],
              id.vars = c("ea","ma","d"),
              variable.name = "Species",
              value.name = "Correlation")

# p2<-ggplot(data=databis)+
#   geom_line(aes(d,Correlation,colour=Species),size=1.5)+
#   corr_colour_TL_4+
#   theme
# legend<-get_legend(p2)

p2<-ggplot(data=databis)+
      geom_line(aes(d,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      theme+#theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle("Weighted dependencies")

# BIOMASS CV
databis<-data_V[,-which(names(data_V)%in%names(params_data)[-which(names(params_data)%in%c("ea","ma","d"))])]
databis<-table_for_matrix(databis,3)
databis<-databis[databis$species_1==databis$species_2 & databis$community_1==databis$community_2,]
databis$species_2<-NULL
databis$community_2<-NULL
names(databis)[names(databis)=="species_1"]="species"
names(databis)[names(databis)=="community_1"]="community"
databis_2<-data_B[,-which(names(data_B)%in%names(params_data)[-which(names(params_data)%in%c("ea","ma","d"))])]
databis_2<-table_for_biomass(databis_2,3)
databis<-merge(databis,databis_2,by=c("ea","ma","d","species","community"))
rm(databis_2)
databis$CV<-sqrt(databis$value)/databis$biomass
databis$species<-as.factor(databis$species)
databis$community<-as.factor(databis$community)

p4<-ggplot(data=databis)+
      geom_line(aes(d,CV,colour=species,linetype=community),size=1.1)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      patch_line+
      theme+
      x_axis_log10_short+
      y_axis_log10_short+
      xlab(label_dispersal)+
      ylab(label_CV)

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(p4, 1, 0, 1, 1)+
  #draw_plot(legend, 2.05, 1.25, 0.05, 0.5)+
  draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 30)
ggsave(paste(path_figure_supp,"supp_env_everything.pdf",sep=""),graph, width = 19, height = 16, device = cairo_pdf)

########################## ----
# REALISTIC CASE - LOCAL CV #### ----
### SIMULATIONS - d - CORRELATION - EVERYTHING - CONTRIBUTION DEPENDING OR NOT ON GROWTH RATE - CORRELATED ENVIRONMENTAL PERTURBATIONS ####
# parameters
nSpecies=4
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(2,1),c(3,1),c(4,1))), # perturbation of predators in patch 1
                         disp=list(c(1,1,1,1)), # dispersal of all species
                         d=d,
                         S0_self=1e3,
                         S0_prey=1e-3,
                         S0_pred=1e3,
                         C0=0,
                         C0_self=1,
                         C0_prey=1,
                         C0_pred=1,
                         correction=correction,
                         model="all")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
params_data$z=1 # environmental perturbation
VE<-matrix(sigma^2,nSpecies*nCommunity,nSpecies*nCommunity) # correlated perturbations
params_data$VE<-list(VE)
# not weighted
params_data$weight=0
# weighted
params_bis<-params_data
params_bis$weight=1
params_data<-rbind(params_data,params_bis)
rm(params_bis)

# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

# BIOMASS CV
databis<-data_V[,-which(names(data_V)%in%names(params_data)[-which(names(params_data)%in%c("ea","ma","d","weight"))])]
databis<-table_for_matrix(databis,4)
databis<-databis[databis$species_1==databis$species_2 & databis$community_1==databis$community_2,]
databis$species_2<-NULL
databis$community_2<-NULL
names(databis)[names(databis)=="species_1"]="species"
names(databis)[names(databis)=="community_1"]="community"
databis_2<-data_B[,-which(names(data_B)%in%names(params_data)[-which(names(params_data)%in%c("ea","ma","d","weight"))])]
databis_2<-table_for_biomass(databis_2,4)
databis<-merge(databis,databis_2,by=c("ea","ma","d","weight","species","community"))
rm(databis_2)
databis$CV<-sqrt(databis$value)/databis$biomass
databis$species<-as.factor(databis$species)
databis$community<-as.factor(databis$community)
databis$weight<-as.factor(databis$weight)
levels(databis$weight)<-c("not weighted","weighted")

p1<-ggplot(data=databis[databis$community=="1",])+
      geom_line(aes(d,CV,colour=species,linetype=weight),size=1.1)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      weight_line+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      y_axis_log10_short+
      xlab(label_dispersal)+
      ylab(label_CV)+
      ggtitle("Patch 1")

p2<-ggplot(data=databis)+
  geom_line(aes(d,CV,colour=species,linetype=weight),size=1.1)+
  corr_colour_TL_4+
  weight_line+
  theme
legend<-get_legend(p2)

p2<-ggplot(data=databis[databis$community=="2",])+
      geom_line(aes(d,CV,colour=species,linetype=weight),size=1.1)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      weight_line+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      y_axis_log10_short+
      xlab(label_dispersal)+
      ylab(label_CV)+
      ggtitle("Patch 2")

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 2.25), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.1, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure_supp,"supp_env_everything_CV_local.pdf",sep=""),graph, width = 18, height = 8, device = cairo_pdf)

########################## ----
# REALISTIC CASE - TOTAL CV #### ----
### SIMULATIONS - d - CORRELATION - PASSIVE DISPERSAL - CORRELATED ENVIRONMENTAL PERTURBATIONS ####
# parameters
nSpecies=4
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(2,1),c(3,1),c(4,1))), # perturbation of predators in patch 1
                         disp=list(c(1,1,1,1)), # dispersal of prey
                         d=d,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=1,
                         C0=1,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="disp_passive")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
params_data$z=1 # environmental perturbation
VE<-matrix(sigma^2,nSpecies*nCommunity,nSpecies*nCommunity) # correlated perturbations
params_data$VE<-list(VE)
params_data<-params_data[params_data$ea==as.character(ea01) & params_data$ma==as.character(ma10),]

# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-table_for_matrix(data_V[,c(which(names(data_V)%in%c("ea","ma","d")),(which(names(data_V)%in%"VE")+1):dim(data_V)[2])],3)
databis<-databis[databis$species_1==databis$species_2,]
databis$species_2<-NULL
databis<-aggregate(databis$value, by=list(ea=databis$ea,
                                          ma=databis$ma,
                                          d=databis$d,
                                          species=databis$species_1), FUN=sum)
names(databis)[names(databis)=="x"]="variance"

databis_2<-table_for_biomass(data_B[,c(which(names(data_B)%in%c("ea","ma","d")),(which(names(data_B)%in%"VE")+1):dim(data_B)[2])],3)
databis_2<-aggregate(databis_2$biomass, by=list(ea=databis_2$ea,
                                                ma=databis_2$ma,
                                                d=databis_2$d,
                                                species=databis_2$species), FUN=sum)
names(databis_2)[names(databis_2)=="x"]="biomass"

databis<-merge(databis,databis_2,by=c("ea","ma","d","species"))
rm(databis_2)
databis$CV<-sqrt(databis$variance)/databis$biomass

p1<-ggplot(data=databis)+
      geom_line(aes(d,CV,colour=species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      xlab(label_dispersal)+
      ylab(label_CV)+
      ggtitle("Passive dispersal")

### SIMULATIONS - d - CORRELATION - EVERYTHING - CONTRIBUTION DEPENDING ON GROWTH RATE - CORRELATED ENVIRONMENTAL PERTURBATIONS ####
# parameters
nSpecies=4
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(2,1),c(3,1),c(4,1))), # perturbation of predators in patch 1
                         disp=list(c(1,1,1,1)), # dispersal of prey
                         d=d,
                         S0_self=1e3,
                         S0_prey=1e-3,
                         S0_pred=1e3,
                         C0=0,
                         C0_self=1,
                         C0_prey=1,
                         C0_pred=1,
                         correction=correction,
                         model="disp_prey")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
params_data$z=1 # environmental perturbation
VE<-matrix(sigma^2,nSpecies*nCommunity,nSpecies*nCommunity) # correlated perturbations
params_data$VE<-list(VE)
params_data$weight=1
params_data<-params_data[params_data$ea==as.character(ea01) & params_data$ma==as.character(ma10),]
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-table_for_matrix(data_V[,c(which(names(data_V)%in%c("ea","ma","d")),(which(names(data_V)%in%"VE")+1):dim(data_V)[2])],3)
databis<-databis[databis$species_1==databis$species_2,]
databis$species_2<-NULL
databis<-aggregate(databis$value, by=list(ea=databis$ea,
                                          ma=databis$ma,
                                          d=databis$d,
                                          species=databis$species_1), FUN=sum)
names(databis)[names(databis)=="x"]="variance"

databis_2<-table_for_biomass(data_B[,c(which(names(data_B)%in%c("ea","ma","d")),(which(names(data_B)%in%"VE")+1):dim(data_B)[2])],3)
databis_2<-aggregate(databis_2$biomass, by=list(ea=databis_2$ea,
                                                ma=databis_2$ma,
                                                d=databis_2$d,
                                                species=databis_2$species), FUN=sum)
names(databis_2)[names(databis_2)=="x"]="biomass"

databis<-merge(databis,databis_2,by=c("ea","ma","d","species"))
rm(databis_2)
databis$CV<-sqrt(databis$variance)/databis$biomass

p2<-ggplot(data=databis)+
      geom_line(aes(d,CV,colour=species),size=1.5)+
      corr_colour_TL_4+
      theme
legend<-get_legend(p2)

p2<-ggplot(data=databis)+
      geom_line(aes(d,CV,colour=species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      xlab(label_dispersal)+
      ylab(label_CV)+
      ggtitle("Weighted dependencies")

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure_supp,"supp_env_everything_CV.pdf",sep=""),graph, width = 14, height = 6, device = cairo_pdf)

########################## ----
### ADDITIONNAL INFORMATION ####
# parameters
nSpecies=4
nCommunity=1
params_data<-expand.grid(pert=list(list(c(1,1),c(2,1),c(3,1),c(4,1))), # perturbation of species 1
                         disp=list(c(0,0,0,0)), # no dispersal
                         d=0,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=1,
                         C0=0,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="all")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
params_data$z=1 # environmental perturbation
VE<-matrix(sigma^2,nSpecies*nCommunity,nSpecies*nCommunity) # correlated perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

# BIOMASS DISTRIBUTION ####
databis<-melt(data_B[,which(names(data_B)%in%c("ea","ma","d","S0_pred","model","B_11",'B_21',"B_31","B_41"))],
              id.vars = c("ea","ma","d","S0_pred","model"),
              variable.name = "species",
              value.name = "biomass")
databis$x<-as.numeric(databis$species)

p1<-ggplot(data=databis)+
      geom_rect(aes(xmin=x-0.5,xmax=x+0.5,ymin=1e-4,ymax=biomass,fill=species))+
      facet_grid(ma~ea, labeller=label_parsed)+
      fill_colour_TL_4+
      y_axis_log10_short+
      coord_flip()+
      theme+
      xlab('Trophic level')+
      ylab("Biomass")

# WEIGHT DEPENDENCY ####
databis<-NULL
for(i in 1:dim(data_B)[1]){
  weight<-get_weight(as.numeric(data_B[i,which(names(data_B)%in%c("B_11","B_21","B_31","B_41"))]),
                     data_B$e[i],data_B$a[i],data_B$m[i],
                     1,1,1,
                     nSpecies)
  weight$e=data_B$e[i]
  weight$a=data_B$a[i]
  weight$m=data_B$m[i]
  databis<-rbind(databis,weight)
}
databis<-merge(databis,data_B[,which(names(data_B)%in%c("e","a","m","ea","ma"))],
               by=c("e","a","m"))

p3<-ggplot(data=databis)+
      geom_col(aes(species,weight,fill=dependency),position="fill")+
      facet_grid(ma~ea, labeller=label_parsed)+
      fill_colour_weight+
      theme+
      xlab('Trophic level')+
      ylab('Weight of dispersal components')

# CORRELATION MATRIX ####
databis<-data_C[,-which(names(data_C)%in%names(params_data)[-which(names(params_data)%in%c("ea","ma"))])]
databis<-table_for_matrix(databis,2)

p4<-ggplot(data=databis)+
      geom_raster(aes(species_2,species_1,fill=value))+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_grad+
      theme_matrix+
      xlab("")+
      ylab("")

# RELATIVE IMPORTANCE OF DISPERSAL ####
nSpecies=4
nCommunity=1
params_data<-expand.grid(d=d,
                         S0_self=1e3,
                         S0_prey=1e-3,
                         S0_pred=1e3,
                         C0=1,
                         C0_self=1,
                         C0_prey=1,
                         C0_pred=1,
                         correction=correction)
params_data<-merge(params_data_original,params_data)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% get_dispersal_importance(params_data,nSpecies,i)
stopCluster(cl)
data<-params_data
data$'1'=0
data$'2'=0
data$'3'=0
data$'4'=0
for (i in 1:dim(params_data)[1]){
  data[i,which(names(data)%in%as.character(1:nSpecies))]<-unlist(results[[i]])
}

databis<-melt(data[,which(names(data)%in%c("ea","ma","d","1",'2',"3","4"))],
              id.vars = c("ea","ma","d"),
              variable.name = "species",
              value.name = "M")

p2<-ggplot(data=databis)+
      geom_line(aes(d,M,colour=species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      x_axis_log10_short+
      theme+
      xlab(label_dispersal)+
      ylab(label_balance)

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(p4, 1, 0, 1, 1)+
  draw_plot_label(c("A","B","C","D"), c(0,0.95,0,0.95), c(2,2,1,1), size = 30)
ggsave(paste(path_figure_supp,"supp_biomass_corr_M_weight.pdf",sep=""),graph, width = 19, height = 16, device = cairo_pdf)

########################## ----
### SIMULATIONS - d - CORRELATION - EVERYTHING - CONTRIBUTION DEPENDING ON GROWTH RATE - CORRELATED ENVIRONMENTAL PERTURBATIONS - DISPERSAL OF SPECIES 4 ####
# parameters
nSpecies=4
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(2,1),c(3,1),c(4,1))), # perturbation of predators in patch 1
                         disp=list(c(0,0,0,1)), # dispersal of prey
                         d=d,
                         S0_self=1e3,
                         S0_prey=1e-3,
                         S0_pred=1e3,
                         C0=0,
                         C0_self=1,
                         C0_prey=1,
                         C0_pred=1,
                         correction=correction,
                         model="disp_prey")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
params_data$z=1 # environmental perturbation
VE<-matrix(sigma^2,nSpecies*nCommunity,nSpecies*nCommunity) # correlated perturbations
params_data$VE<-list(VE)
params_data$weight=1
params_data<-params_data[params_data$ea==as.character(ea01) & params_data$ma==as.character(ma10),]
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","C_11_12",'C_21_22','C_31_32','C_41_42'))],
              id.vars = c("ea","ma","d"),
              variable.name = "species",
              value.name = "correlation")

p1<-ggplot(data=databis)+
      geom_line(aes(d,correlation,colour=species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle("Dispersal of species 4")

### SIMULATIONS - d - CORRELATION - EVERYTHING - CONTRIBUTION DEPENDING ON GROWTH RATE - CORRELATED ENVIRONMENTAL PERTURBATIONS - DISPERSAL OF SPECIES 4 AND 3 ####
# parameters
nSpecies=4
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(2,1),c(3,1),c(4,1))), # perturbation of predators in patch 1
                         disp=list(c(0,0,1,1)), # dispersal of prey
                         d=d,
                         S0_self=1e3,
                         S0_prey=1e-3,
                         S0_pred=1e3,
                         C0=0,
                         C0_self=1,
                         C0_prey=1,
                         C0_pred=1,
                         correction=correction,
                         model="disp_prey")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
params_data$z=1 # environmental perturbation
VE<-matrix(sigma^2,nSpecies*nCommunity,nSpecies*nCommunity) # correlated perturbations
params_data$VE<-list(VE)
params_data$weight=1
params_data<-params_data[params_data$ea==as.character(ea01) & params_data$ma==as.character(ma10),]
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","C_11_12",'C_21_22','C_31_32','C_41_42'))],
              id.vars = c("ea","ma","d"),
              variable.name = "species",
              value.name = "correlation")

p2<-ggplot(data=databis)+
      geom_line(aes(d,correlation,colour=species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle("Dispersal of species 4 and 3")

### SIMULATIONS - d - CORRELATION - EVERYTHING - CONTRIBUTION DEPENDING ON GROWTH RATE - CORRELATED ENVIRONMENTAL PERTURBATIONS - DISPERSAL OF SPECIES 4, 3 AND 2 ####
# parameters
nSpecies=4
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1),c(2,1),c(3,1),c(4,1))), # perturbation of predators in patch 1
                         disp=list(c(0,1,1,1)), # dispersal of prey
                         d=d,
                         S0_self=1e3,
                         S0_prey=1e-3,
                         S0_pred=1e3,
                         C0=0,
                         C0_self=1,
                         C0_prey=1,
                         C0_pred=1,
                         correction=correction,
                         model="disp_prey")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
params_data$z=1 # environmental perturbation
VE<-matrix(sigma^2,nSpecies*nCommunity,nSpecies*nCommunity) # correlated perturbations
params_data$VE<-list(VE)
params_data$weight=1
params_data<-params_data[params_data$ea==as.character(ea01) & params_data$ma==as.character(ma10),]
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","C_11_12",'C_21_22','C_31_32','C_41_42'))],
              id.vars = c("ea","ma","d"),
              variable.name = "species",
              value.name = "correlation")

p3<-ggplot(data=databis)+
      geom_line(aes(d,correlation,colour=species),size=1.5)+
      corr_colour_TL_4+
      theme
legend<-get_legend(p3)

p3<-ggplot(data=databis)+
      geom_line(aes(d,correlation,colour=species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_4+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle("Dispersal of species 4, 3 and 2")

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 3.15), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(p3, 2, 0, 1, 1)+
  draw_plot(legend, 3.05, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B","C"), c(0,1,2), c(1,1,1), size = 30)
ggsave(paste(path_figure_supp,"supp_env_everything_disp.pdf",sep=""),graph, width = 18, height = 6, device = cairo_pdf)

########################## ----
# LINEAR DENSITY DEPENDENCY #### ----
# DISPERSAL DEPENDING ON SELF INTERACTIONS - CORRELATION ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         d=d,
                         S0_self=0,
                         S0_prey=0,
                         S0_pred=0,
                         C0=0,
                         C0_self=1,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="disp_pred")
params_data<-rbind(params_data,expand.grid(pert=list(list(c(2,1))), # perturbation of predators in patch 1
                                           disp=list(c(1,0)), # dispersal of prey
                                           d=d,
                                           S0_self=c(1e-3,1e3),
                                           S0_prey=1,
                                           S0_pred=1,
                                           C0=0,
                                           C0_self=1,
                                           C0_prey=0,
                                           C0_pred=0,
                                           correction=correction,
                                           model="disp_prey"))
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_linear,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","model","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","model"),
              variable.name = "Species",
              value.name = "Correlation")

p1<-ggplot(data=databis[databis$model=="disp_pred",])+
      geom_line(aes(d,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle(expression(atop("Dispersal of predators","depending on their density")))

# DISPERSAL OF PREY DEPENDING ON PRED - d - CORRELATION ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(2,1))), # perturbation of predators in patch 1
                         disp=list(c(1,0)), # dispersal of prey
                         d=d,
                         S0_self=0,
                         S0_prey=0,
                         S0_pred=0,
                         C0=0,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=1,
                         correction=correction,
                         model="disp_prey")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_linear,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d"),
              variable.name = "Species",
              value.name = "Correlation")

p2<-ggplot(data=databis)+
      geom_line(aes(d,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle(expression(atop("Dispersal of prey","depending on predators")))

# DISPERSAL OF PRED DEPENDING ON PREY - d - CORRELATION ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         d=d,
                         S0_self=0,
                         S0_prey=0,
                         S0_pred=0,
                         C0=2,
                         C0_self=0,
                         C0_prey=1,
                         C0_pred=0,
                         correction=correction,
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_linear,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","C0","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","C0"),
              variable.name = "Species",
              value.name = "Correlation")

p3<-ggplot(data=databis)+
      geom_vline(xintercept=1e3,linetype='dashed')+
      geom_line(aes(d,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(label_dispersal)+
      ylab(label_correlation)+
      ggtitle(expression(atop("Dispersal of predators","depending on prey ("*italic(C["0"])*"=2)")))

# DISPERSAL OF PRED DEPENDING ON PREY - C0 - CORRELATION ####
# parameters
nSpecies=2
nCommunity=2
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of prey in patch 1
                         disp=list(c(0,1)), # dispersal of predators
                         d=1e3,
                         S0_self=0,
                         S0_prey=0,
                         S0_pred=0,
                         C0=seq(0.5,110,0.1),
                         C0_self=0,
                         C0_prey=1,
                         C0_pred=0,
                         correction=correction,
                         model="disp_pred")
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_linear,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-melt(data_C[,which(names(data_C)%in%c("ea","ma","d","C0","C_11_12",'C_21_22'))],
              id.vars = c("ea","ma","d","C0"),
              variable.name = "Species",
              value.name = "Correlation")

p4<-ggplot(data=databis)+
      geom_line(aes(C0,Correlation,colour=Species),size=1.5)+
      corr_colour_TL_2+
      theme
legend<-get_legend(p4)

p4<-ggplot(data=databis)+
      geom_vline(xintercept=2,linetype='dashed')+
      geom_line(aes(C0,Correlation,colour=Species),size=1.5)+
      facet_grid(ma~ea, labeller=label_parsed)+
      corr_colour_TL_2+
      theme+theme(legend.position = "none")+
      x_axis_log10_short+
      ylim(-1,1)+
      xlab(expression("Basal dispersal intensity "*italic(C["0"])))+
      ylab(label_correlation)+
      ggtitle(expression(atop("Dispersal of predators","depending on prey ("*italic(d["i"])*"=10"^"3"*")")))

### FINAL GRAPH ####
graph<-ggdraw(xlim = c(0, 2.15), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p2, 1, 1, 1, 1)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot(p4, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 1.25, 0.05, 0.5)+
  draw_plot(legend, 2.05, 0.25, 0.05, 0.5)+
  draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 30)
ggsave(paste(path_figure_supp,"supp_linear.pdf",sep=""),graph, width = 15, height = 16, device = cairo_pdf)

########################## ----
# VARIANCE TRANSMISSION #### ----
# parameters - 4 species
nSpecies=4
nCommunity=1
params_data<-expand.grid(pert=list(list(c(1,1))), # perturbation of species 1 in patch 1
                         disp=list(c(0,0,0,0)), # no dispersal
                         d=0,
                         S0_self=1,
                         S0_prey=1,
                         S0_pred=1,
                         C0=0,
                         C0_self=0,
                         C0_prey=0,
                         C0_pred=0,
                         correction=correction,
                         model="1")
params_data<-rbind(params_data,expand.grid(pert=list(list(c(4,1))), # perturbation of species 4 in patch 1
                                           disp=list(c(0,0,0,0)), # no dispersal
                                           d=0,
                                           S0_self=1,
                                           S0_prey=1,
                                           S0_pred=1,
                                           C0=0,
                                           C0_self=0,
                                           C0_prey=0,
                                           C0_pred=0,
                                           correction=correction,
                                           model="4"))
params_data<-rbind(params_data,expand.grid(pert=list(list(c(1,1),c(2,1),c(3,1),c(4,1))), # perturbation of species 4 in patch 1
                                           disp=list(c(0,0,0,0)), # no dispersal
                                           d=0,
                                           S0_self=1,
                                           S0_prey=1,
                                           S0_pred=1,
                                           C0=0,
                                           C0_self=0,
                                           C0_prey=0,
                                           C0_pred=0,
                                           correction=correction,
                                           model="1234"))
params_data<-merge(params_data_original,params_data)
for(i in 1:dim(params_data)[1]){
  params_data$disp[[i]]=params_data$disp[[i]]*params_data$d[i]
}
VE<-diag(rep(sigma^2,nSpecies*nCommunity)) # independent perturbations
params_data$VE<-list(VE)
# simulations - 4 species
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  
results<-foreach(i=1:dim(params_data)[1]) %dopar% analytic(params_data,jacobian_disp_hauzy,nSpecies,nCommunity,i)
stopCluster(cl)
list_data<-create_data(params_data,nSpecies,nCommunity)
data_B<-list_data[[1]][,-which(names(list_data[[1]])%in%c("pert","disp"))]
data_V<-list_data[[2]][,-which(names(list_data[[2]])%in%c("pert","disp"))]
data_C<-list_data[[3]][,-which(names(list_data[[3]])%in%c("pert","disp"))]
for (i in 1:dim(params_data)[1]){
  data_B[i,which(names(data_B)%in%list_data[[4]])]<-unlist(results[[i]][1])
  data_V[i,which(names(data_V)%in%list_data[[5]])]<-unlist(results[[i]][2])
  data_C[i,which(names(data_C)%in%list_data[[6]])]<-unlist(results[[i]][3])
}

databis<-data_V[,-which(names(data_V)%in%names(params_data)[-which(names(params_data)%in%c("ea","ma","sigma","model"))])]
databis<-table_for_matrix(databis,4)
databis<-databis[databis$species_1==databis$species_2 & databis$community_1==databis$community_2,]
databis$species_2<-NULL
databis$community_2<-NULL
names(databis)[names(databis)=="species_1"]="species"
names(databis)[names(databis)=="community_1"]="community"
names(databis)[names(databis)=="value"]="variance"
databis$x<-as.numeric(databis$species)

p1<-ggplot(data=databis[databis$model=="1",])+
  geom_rect(aes(xmin=x-0.5,xmax=x+0.5,ymin=1e-16,ymax=variance,fill=species))+
  facet_grid(ma~ea, labeller=label_parsed)+
  fill_colour_TL_4+
  y_axis_log10_short+
  theme+theme(legend.position = "none")+
  xlab('Trophic level')+
  ylab("Variance")+
  ggtitle("Species 1 perturbed")

p2<-ggplot(data=databis[databis$model=="4",])+
  geom_rect(aes(xmin=x-0.5,xmax=x+0.5,ymin=1e-16,ymax=variance,fill=species))+
  facet_grid(ma~ea, labeller=label_parsed)+
  fill_colour_TL_4+
  y_axis_log10_short+
  theme+theme(legend.position = "none")+
  xlab('Trophic level')+
  ylab("Variance")+
  ggtitle("Species 4 perturbed")

p3<-ggplot(data=databis[databis$model=="1234",])+
  geom_rect(aes(xmin=x-0.5,xmax=x+0.5,ymin=1e-16,ymax=variance,fill=species))+
  facet_grid(ma~ea, labeller=label_parsed)+
  fill_colour_TL_4+
  y_axis_log10_short+
  theme+theme(legend.position = "none")+
  xlab('Trophic level')+
  ylab("Variance")+
  ggtitle("All species perturbed")

databis<-params_data_original
databis$x=0.5
databis$y=0.5

p4<-ggplot(data=databis)+
  geom_raster(aes(x,y,fill=m))+
  facet_grid(ma~ea, labeller=label_parsed)+
  scale_fill_viridis(name=expression(paste("Metabolic\nrate ratio",italic(m))),
                     breaks = unique(databis$m),
                     trans = "log10")+
  theme_matrix+theme(axis.text=element_blank(),
                     legend.key=element_blank(),
                     legend.text.align = 1)+
  scale_x_continuous(trans=log10_trans())+
  scale_y_continuous(trans=log10_trans())+
  xlab("")+
  ylab("")