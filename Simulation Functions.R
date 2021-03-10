

simpo = function(sp_coef, quadsenv, intsamp=-1.5, npoint_samp=150,
                  samp=c(TRUE, FALSE), plot=c(TRUE, FALSE), 
                  compPlot=c(TRUE, FALSE)){
  
  # Use the env_mat to create window and quadrature
  #winMix   = owin(xrange = c(min(quadsenv$X)-0.5, max(quadsenv$X)+0.5), 
  #                yrange = c(min(quadsenv$Y)-0.5, max(quadsenv$Y)+0.5))
  #quads. = ppp(quadsenv$X, quadsenv$Y, window = winMix)
  #Qmask = makeMask(quads.)
  
  # prepare the environmental matrix with intercept column
  env_mat = as.matrix(cbind(rep(1,dim(quadsenv)[1]), quadsenv[,3:dim(quadsenv)[2]]))
  
  # Species true intensity
  # env_mat needs to have X and Y too
  sp_int = exp(env_mat[,-c(dim(env_mat)[2])] %*% sp_coef[-c(dim(env_mat)[2])])
  
  LVint = levelplot(sp_int ~ quadsenv$X +quadsenv$Y, asp = "iso", main = "True intensity",
                    col.regions=viridis(30),ylab="", xlab="")
  
  sp_int_im = as.im(data.frame(x = quadsenv$X, y = quadsenv$Y, z = sp_int))
  
  # Generate points from the true intensity
  species_sim = rpoispp(sp_int_im)

  quads = data.frame(quadsenv)
  
  species_xy = data.frame(X = species_sim$x, Y = species_sim$y)
  
  species_env = newenv.var(sp.xy = species_xy, env.grid = quads, 
                      env.scale = 0.5, coord = c("X", "Y"), 
                      file.name = NA)
  
  if(isFALSE(samp)){
    if(isTRUE(plot)){
      xy1 = xyplot(species_sim$y~species_sim$x, ylab="", xlab="", col="navy")
      
      comb_levObj <- c(LVint, xy1, layout = c(2, 1), merge.legends = T)
      update(comb_levObj, main="Species intensity and true species pattern")
    }
  }else{
    po_int = exp(env_mat %*% sp_coef)
    LVpo = levelplot(po_int ~ quadsenv$X + quadsenv$Y, asp = "iso", main = "sampled intensity",
              col.regions=viridis(30),ylab="", xlab="")
    
    
    # point pattern sampling according to the covariates
    podata_X = data.frame(Intercept = 1, species_env[,c(3:dim(quads)[2])])
    podata_beta = c(intsamp, sp_coef[-1])
    podata_intensity = exp(as.matrix(podata_X) %*% podata_beta)
    
    POdata_rows = sample(1:species_sim$n, npoint_samp, prob = podata_intensity)
    POdata = species_sim[POdata_rows]

    
    xyplot(POdata$y~POdata$x, ylab="", xlab="", pch=16,
           col=rgb(1, 0.5, 0, alpha=0.4), cex=1.2)
    
    if(isTRUE(compPlot)){
      xy = xyplot(POdata$y~POdata$x, ylab="", xlab="", pch=16,
             col=rgb(1, 0.5, 0, alpha=0.4), cex=1.2)
      
      comb_levObj <- c(LVpo , xy, layout = c(2, 1), merge.legends = T)
      update(comb_levObj, main="Species intensity and sampled pattern")
    }
  }
  
  title = c("Species intensity and true species pattern", "Species intensity and sampled pattern")
  
  if(isTRUE(compPlot)){
    if(isTRUE(samp)){
      i=2
      return(list(sp=species_sim, sp_int=sp_int, PO=POdata, po_intensity=podata_intensity,
                  sp_coef=sp_coef, quadsenv=quadsenv,
                  update(comb_levObj, main=title[i])))
    }else{
      i=1
      return(list(sp=species_sim, sp_int=sp_int,
                  sp_coef=sp_coef, quadsenv=quadsenv,
                  update(comb_levObj, main=title[i])))
    }
  }else{
    comb_levObj = NA
  }
  
  
}


simocc=function(sp_coef, quadsenv, upres, vardetect, n_visits=5,
                   occ_plot=c(TRUE, FALSE)){
  # prepare the environmental matrix with intercept column
  env_mat = as.matrix(cbind(rep(1,dim(quadsenv)[1]), quadsenv[,3:dim(quadsenv)[2]]))
  
  # Species true intensity
  # env_mat needs to have X and Y too
  sp_int = exp(env_mat %*% sp_coef)
  
  LVint = levelplot(sp_int ~ quadsenv$X +quadsenv$Y, asp = "iso", main = "True intensity",
                    col.regions=viridis(30),ylab="", xlab="")
  
  sp_int_im = as.im(data.frame(x = quadsenv$X, y = quadsenv$Y, z = sp_int))
  
  # Generate points from the true intensity
  sp_sim = rpoispp(sp_int_im)
  xyplot(sp_sim$y~sp_sim$x, ylab="", xlab="", col="navy")
  
  # *Generate occupancy data from species true intensity*
  # Creating the grid sites
  XY.occ = expand.grid(seq(0, max(quadsenv$X), max(quadsenv$X)*10/100), 
                       seq(0, max(quadsenv$Y), max(quadsenv$Y)*10/100)) 
  X.occ = XY.occ[,1]
  Y.occ = XY.occ[,2]
  
  # Use the env_mat to create window and quadrature
  winsim  = owin(xrange = c(min(quadsenv$X)-0.5, max(quadsenv$X)+0.5), 
                  yrange = c(min(quadsenv$Y)-0.5, max(quadsenv$Y)+0.5))
  
  sp_and_occ = ppp(x = c(X.occ, sp_sim$x), y = c(Y.occ, sp_sim$y),
                   marks = c(rep("Occ", length(X.occ)), 
                             rep("Sp", sp_sim$n)),
                   window = winsim)
  dist_sp_occ = nndist(sp_and_occ, k = 1, by = as.factor(marks(sp_and_occ)))
  
  # compute distance to nearest species
  sp_occ_dists = dist_sp_occ[1:length(X.occ), 2]
  
  # species considered present at site if nearest one is within 0.25 units
  occ_present = as.numeric(sp_occ_dists <= upres)
  table(occ_present)
  Tocc = table(occ_present) 
  occ_present2 = occ_present
  
  occ_paste = paste(X.occ, Y.occ)
  quad_paste = paste(quad$X, quad$Y)
  occ_quadrow = match(occ_paste, quad_paste)
  
  
  vdetect = quadsenv[,vardetect]
  
  p_detect = clogloginv(vdetect[occ_quadrow])*occ_present  
  
  sim_history = matrix(as.integer(matrix(rep(p_detect, times = n_visits),
                                         length(X.occ), n_visits) >
                                    matrix(runif(length(X.occ)*n_visits),
                                           length(X.occ), n_visits)),
                       length(X.occ), n_visits)
  
  # if(isTRUE(v_plot)){
  #   # Plot of 4 out the 5 visits from the history
  #   
  #   L= list()
  #   for (i in 1:n_visits) {
  #     L[[i]] = sim_history[,i]
  #   }
  #   
  #   Grplot = Groupplot(listLplot=L, X.occ, Y.occ, Ncomp=ceiling(length(L)))
  #   
  #  
  # }
  
  if(isTRUE(occ_plot)){
    # Occupancy plot
    site_sum = apply(sim_history, 1, sum)
    plot_xy = expand.grid(seq(0, max(quadsenv$X), (max(quadsenv$X)*10/100)/2), 
                         seq(0, max(quadsenv$Y), (max(quadsenv$Y)*10/100)/2))
    plot_x = plot_xy[,1]
    plot_y = plot_xy[,2]
    plot_id = paste(plot_x, plot_y)
    occ_id = paste(X.occ, Y.occ)
    occ_match = match(occ_id, plot_id)
    plot_z = rep(NA, length(plot_x))
    plot_z[occ_match] = site_sum
    Lz = levelplot(plot_z ~ plot_x + plot_y, asp = "iso", cuts = 5,
                   col.regions = rainbow(6, start = 0.03, end = 0.17)[6:1])
  }
  
  
  if(isFALSE(occ_plot)){
    return(list(Tocc, sp_sim, sp_coef=sp_coef, quadsenv=quadsenv))
  }
  if(isTRUE(occ_plot)){
    return(list(Tocc, sp_sim, sp_coef=sp_coef, quadsenv=quadsenv, Lz))
  }
  #if(isTRUE(v_plot)& isTRUE(occ_plot)){
  #  return(list(Tocc, sp_sim, Grplot, sp_coef=sp_coef, quadsenv=quadsenv, Lz))
  #}
}



misidpoints = function(simpoRes, perc_mis, compPlot=c(TRUE, FALSE)){
  
  POdata = simpoRes$PO
  quads = data.frame(simpoRes$quadsenv)
  
  
  xy = xyplot(POdata$y~POdata$x, pch=16, cex=1.2,
                col=rgb(0.6, 0.8, 0.3, alpha=0.4))
  
  # Misidentification
      # Select a PO pattern with false negative
      spPOb.xy = data.frame(X=POdata$x, Y=POdata$y)
      
      sp_envb = newenv.var(sp.xy = spPOb.xy, env.grid = quads, 
                           env.scale = 0.5, coord = c("X", "Y"), file.name = NA)
      pomis_X = data.frame(Intercept = 1, sp_envb[,c(3:dim(quads)[2])])
      po_beta = simpoRes$sp_coef
      pomis_intensity = exp(as.matrix(pomis_X) %*% po_beta)
      
      POmis_rows = sample(1:POdata$n,  perc_mis*POdata$n, prob = pomis_intensity)
      POmis = POdata[POmis_rows]
      
      posup = superimpose(PO=POdata, POmis=POmis)
      
      xmis = xyplot(posup$y~posup$x, pch=16,
             groups = posup$marks, cex=1.2,
             col=c(rgb(0.6, 0.8, 0.3, alpha=0.4),rgb(0.3, 0.7, 0.8, alpha=0.4)))
      
      comb_levObj <- c(xy, xmis, layout = c(2, 1), merge.legends = T)

  
  title = c("species sampled and misidentification")
  
  if(isTRUE(compPlot)){
      return(list(PO=POdata, 
                  POmis = POmis, pomis_intensity=pomis_intensity, 
                  sp_coef=po_beta, quadsenv=quads,
                  update(comb_levObj, main=title)))
  }else{
    comb_levObj = NA
  }
  
  
}


multiplepo = function(sp_coefmult, quadsenv, intsampmult, npoint_sampmult,
                     plot=c(TRUE, FALSE), 
                     compPlot=c(TRUE, FALSE)){
  
  # prepare the environmental matrix with intercept column
  env_mat = as.matrix(cbind(rep(1,dim(quadsenv)[1]), quadsenv[,3:dim(quadsenv)[2]]))
  
  
  sp_int = LVint = sp_int_im = sp_sim = sp_xy = sp_env = 
    po_beta = po_intensity= po_X = PO_rows =PO = xy1 =xy = 
    coordx.list = coordy.list = mark.list = list()
  for (i in 1:length(intsampmult)) {
    # Species true intensity
    # env_mat needs to have X and Y too
    sp_int[[i]] = exp(env_mat %*% sp_coefmult[[i]])
    
    LVint[[i]] = levelplot(sp_int[[i]] ~ quadsenv$X +quadsenv$Y, asp = "iso", main = "True intensity",
                      col.regions=viridis(30),ylab="", xlab="")
    
    sp_int_im[[i]] = as.im(data.frame(x = quadsenv$X, y = quadsenv$Y, z = sp_int[[i]]))
    
    # Generate points from the true intensity
    sp_sim[[i]] = rpoispp(sp_int_im[[i]])
    
    quads = data.frame(quadsenv)
    
    sp_xy[[i]] = data.frame(X = sp_sim[[i]]$x, Y = sp_sim[[i]]$y)
    
    sp_env[[i]] = newenv.var(sp.xy = sp_xy[[i]], env.grid = quads, 
                        env.scale = 0.5, coord = c("X", "Y"), 
                        file.name = NA)
    
    
    # point pattern sampling according to the covariates
    po_X[[i]] = data.frame(Intercept = 1, sp_env[[i]][,c(3:dim(quadsenv)[2])])
    po_beta[[i]] = c(intsampmult[[i]], sp_coefmult[[i]][-1])
    po_intensity[[i]] = exp(as.matrix(po_X[[i]]) %*% po_beta[[i]])
    PO_rows[[i]] = sample(1:sp_sim[[i]]$n, npoint_sampmult[[i]], prob = po_intensity[[i]])
    PO[[i]] = sp_sim[[i]][PO_rows[[i]]]
    
    coordx.list[[i]] = PO[[i]]$x
    coordy.list[[i]] = PO[[i]]$y
    mark.list[[i]] = rep(paste("Sp", i, sep = ""), PO[[i]]$n)

  }
  colvec = c(rgb(1,0.6,0, 0.4), rgb(0.2,0.4,0.6, 0.4), rgb(0.6,0.2,0.8, 0.4), rgb(1,0,0.6, 0.4), 
             rgb(0.2,0.8,0.2, 0.4), rgb(0.4,0,0, 0.4))

  winsim  = owin(xrange = c(min(quadsenv$X)-0.5, max(quadsenv$X)+0.5), 
                 yrange = c(min(quadsenv$Y)-0.5, max(quadsenv$Y)+0.5))
  
  all_po = ppp(x = c(unlist(coordx.list)), 
                 y = c(unlist(coordy.list)), window = winsim,
                 marks = as.factor(c(unlist(mark.list))))
  
  marksvec = c("po1", "po2", "po3", "po4", "po5", "po6")
  
  xy = xyplot(all_po$y~all_po$x, ylab="", xlab="", pch=16,
         col=colvec[1:length(table(all_po$marks))], cex=1.2,
         groups=all_po$marks, auto.key=list(space="right",text=marksvec[1:length(table(all_po$marks))],
                                            cex=1.8,points=FALSE,
                                            col=colvec[1:length(table(all_po$marks))]))

  return(list(sp=sp_sim, sp_int=sp_int, PO=PO, po_intensity=po_intensity,
              all_po = all_po, xy = xy))

  
}


deathpoints = function(simpoRes, prob_death, sd=0.01, random=FALSE){
  
  if(isFALSE(random)){
    PO_xy = data.frame(X = simpoRes$PO$x, Y = simpoRes$PO$y)
    
    quads = data.frame(simpoRes$quadsenv)
    PO_env = newenv.var(sp.xy = PO_xy, env.grid = quads, 
                        env.scale = 0.5, coord = c("X", "Y"), file.name = NA)
    
    po_X = data.frame(Intercept = 1, PO_env[,c(3:dim(quads)[2])])
    po_beta = simpoRes$sp_coef
    po_intensity = exp(as.matrix(po_X) %*% po_beta)
    
    # Species intensity is used to consider the more suitable habitat
    P_retention = as.matrix(po_intensity)
    notdead = sample(1:simpoRes$PO$n, round(rnorm(1, 1-prob_death, sd)*simpoRes$PO$n),
                     prob = P_retention)
    PO_dead = simpoRes$PO[-notdead]
    PO_dead
    
    plot.new()
    plot(simpoRes$PO, col="white", main= "Death process: habitat suitability")
    points(simpoRes$PO, pch=16, col=rgb(0,0,0.6, 0.4), cex=1.8)
    points(PO_dead, pch=16, col=rgb(0,1,0.8, 0.4), cex=1.8)
    
    return(list(PO=simpoRes$PO, PO_dead=PO_dead, prob_death=prob_death))
  }else{
    
    died = sample(1:simpoRes$PO$n, round(rnorm(1, prob_death, sd)*simpoRes$PO$n))
    PO_yr2 = simpoRes$PO[-died]
    PO_died = simpoRes$PO[died]
    
    plot(simpoRes$PO, col="white", main="Death process: random")
    points(PO_yr2, pch=16, col=rgb(0,0,0.6, 0.4), cex=1.8)
    points(PO_died, pch=16, col=rgb(0,1,0.8, 0.4), cex=1.8)
  }
   
}


birthpoints=function(simpoRes, prop_offsp, radius=1, birth=c(TRUE,FALSE)){
  
  PO_xy = data.frame(X = simpoRes$PO$x, Y = simpoRes$PO$y)
  
  quads = data.frame(simpoRes$quadsenv)
  PO_env = newenv.var(sp.xy = PO_xy, env.grid = quads, 
                      env.scale = 0.5, coord = c("X", "Y"), file.name = NA)
  
  po_X = data.frame(Intercept = 1, PO_env[,c(3:dim(quads)[2])])
  po_beta = simpoRes$sp_coef
  po_intensity = exp(as.matrix(po_X) %*% po_beta)
  
  # Species intensity is used to consider the more suitable habitat
  P_retention = as.matrix(po_intensity)
  
  
  if(isTRUE(birth)){
    Parents_rows = c() #vector indicating the sampled rows
    Xparents_add = sample(1:simpoRes$PO$n, prop_offsp*simpoRes$PO$n, prob = P_retention)
    Parents_rows = c(Parents_rows, Xparents_add)
    Parents = simpoRes$PO[Parents_rows]
    
    Cparents = cbind(Parents$x, Parents$y)
    
    pt.disc = coord.pt = coordallx = coordally = list()
    for (i in 1:Parents$n) {
      pt.disc[[i]] = runifdisc(n=1, radius, centre=Cparents[i,], mask=T)
      
      coord.pt[[i]] = as.data.frame(cbind(pt.disc[[i]]$x, pt.disc[[i]]$y))
      colnames(coord.pt[[i]])=c("x", "y")
      coordallx[[i]] = coord.pt[[i]]$x
      coordally[[i]] = coord.pt[[i]]$y	
    }
    
    winsim  = owin(xrange = c(min(simpoRes$quadsenv$X)-0.5, max(simpoRes$quadsenv$X)+0.5), 
                   yrange = c(min(simpoRes$quadsenv$Y)-0.5, max(simpoRes$quadsenv$Y)+0.5))
    
    Offsprg = ppp(x = c(unlist(coordallx)), 
                  y = c(unlist(coordally)), window = winsim)
    
    # We gather the original and offspring patterns to get the PO data for the year 2.
    PO_year2 = superimpose(simpoRes$PO, Offsprg)
    
    plot(PO_year2, col="white", main="PO with offsprings Year2")
    points(simpoRes$PO, col=rgb(0.4,0,0.6, 0.4), pch=16, cex=1.8)
    points(Offsprg, col=rgb(1,0.6,0.7, 0.6), pch=15, cex=1.8)
    
    return(list(PO=simpoRes$PO, po_intensity=simpoRes$po_intensity, Offsprg =Offsprg, 
                PO_year2=PO_year2, prop_offsp=prop_offsp))
    
  }else{
    PO_rows = sample(1:simpoRes$PO$n, prop_offsp*simpoRes$PO$n, prob = P_retention)
    newies = simpoRes$PO[PO_rows]
    
    POb_yr2 = simpoRes$PO
    POb_yr2$x = c(POb_yr2$x, newies$x)
    POb_yr2$y = c(POb_yr2$y, newies$y)
    
    plot(simpoRes$PO, col="white", main ="PO and colonized points")
    points(simpoRes$PO[-PO_rows], pch=16, col=rgb(0.4,0,0.6, 0.4), cex=1.8)
    points(newies, pch=15, col=rgb(1,0.6,0.7, 0.6), cex=1.8)
    
    return(list(PO=simpoRes$PO, po_intensity=simpoRes$po_intensity, newies =newies, 
                POb_yr2=POb_yr2, prop_offsp=prop_offsp))
  }

}


movepoints=function(simpoRes, mov_type = c("random", "intensity", "center"),p_scale= NULL,
                 Allpt = c(TRUE, FALSE), mean=0, sd=2, n_var=NULL, p_move=NULL){
  
  mov_type = match.arg(mov_type)
  
  if(mov_type == "random"){
    # Case 1: movement from randorm generation using the normal function
    
    if(isTRUE(Allpt)){
      sp1_mov_yr2 = simpoRes$PO
      sp1_mov_yr2$x = simpoRes$PO$x + rnorm(simpoRes$PO$n, mean, sd)
      sp1_mov_yr2$y = simpoRes$PO$y + rnorm(simpoRes$PO$n, mean, sd)
      
      par(mfrow=c(1,2))
      plot(simpoRes$PO, main="before", col="white")
      points(simpoRes$PO, pch=16, col=rgb(0.8,0,0, 0.4))
      
      plot(sp1_mov_yr2, main="after", pch=16, col="white")
      points(sp1_mov_yr2, pch=16, col=rgb(0,0.8,1, 0.4))
      
      par(mfrow=c(1,1))
    }else{
      # choose which points to move
      sp1_sub = sample(1:simpoRes$PO$n, floor(p_move*simpoRes$PO$n))
      sp1_sub_static  = simpoRes$PO[-sp1_sub]
      
      sp1_movsub_yr2 = simpoRes$PO[sp1_sub]
      sp1_movsub_yr2b = simpoRes$PO[sp1_sub]

      sp1_movsub_yr2$x = sp1_movsub_yr2$x + rnorm(sp1_movsub_yr2$n, mean, sd)
      sp1_movsub_yr2$y = sp1_movsub_yr2$y + rnorm(sp1_movsub_yr2$n, mean, sd)
      
      sp1_mov_yr2= superimpose(static=sp1_sub_static, move=sp1_movsub_yr2)
      
      par(mfrow=c(1,2))
      plot(simpoRes$PO, main="before", col="white")
      points(simpoRes$PO, pch=16, col=rgb(0.8,0,0, 0.4))
      points(sp1_movsub_yr2b, pch=16, col=rgb(0,0.8,1, 0.4))
      
      plot(sp1_mov_yr2, main="after", pch=c(16,16), cols=c(rgb(0.8,0,0, 0.4), rgb(0,0.8,1, 0.4)), 
           legend=F)
      
      par(mfrow=c(1,1))
    }
    
    return(list(PO=simpoRes$PO, PO_mov=sp1_mov_yr2))
  }
    
  if(mov_type == "intensity"){
    # Case 2: movement in more suitable environment
    if(isTRUE(Allpt)){
      sp1_mov_yr2 = simpoRes$PO
      sp1_mov_yr2$x = simpoRes$PO$x + rnorm(simpoRes$PO$n, mean, 
                                                max(simpoRes$po_intensity)/range(simpoRes$po_intensity))
      sp1_mov_yr2$y = simpoRes$PO$y + rnorm(simpoRes$PO$n, mean, 
                                                max(simpoRes$po_intensity)/range(simpoRes$po_intensity))
      
      par(mfrow=c(1,2))
      plot(simpoRes$PO, main="before", col="white")
      points(simpoRes$PO, pch=16, col=rgb(0.8,0,0, 0.4))
      
      plot(sp1_mov_yr2, main="after", col="white")
      points(sp1_mov_yr2, pch=16, col=rgb(0,0.8,1, 0.4))
      
      par(mfrow=c(1,1))
    }else{
      sp1_sub = sample(1:simpoRes$PO$n, floor(p_move*simpoRes$PO$n))
      sp1_sub_static  = simpoRes$PO[-sp1_sub]
      
      sp1_movsub_yr2 = simpoRes$PO[sp1_sub]
      sp1_movsub_yr2b = simpoRes$PO[sp1_sub]
      sp1_movsub_yr2$x = sp1_movsub_yr2$x + rnorm(sp1_movsub_yr2$n, mean, 
                                                max(simpoRes$po_intensity)/range(simpoRes$po_intensity))
      sp1_movsub_yr2$y = sp1_movsub_yr2$y + rnorm(sp1_movsub_yr2$n, mean, 
                                                max(simpoRes$po_intensity)/range(simpoRes$po_intensity))
    
      sp1_mov_yr2= superimpose(static=sp1_sub_static, move=sp1_movsub_yr2)
      
      par(mfrow=c(1,2))
      plot(simpoRes$PO, main="before", col="white")
      points(simpoRes$PO, pch=16, col=rgb(0.8,0,0, 0.4))
      points(sp1_movsub_yr2b, pch=16, col=rgb(0,0.8,1, 0.4))
      
      plot(sp1_mov_yr2, main="after", pch=c(16,16), cols=c(rgb(0.8,0,0, 0.4), rgb(0,0.8,1, 0.4)), 
           legend=F)
      
      par(mfrow=c(1,1))
    }
    
    
    return(list(PO=simpoRes$PO, PO_mov=sp1_mov_yr2))
  }
  
  if(mov_type == "center"){
    # Case 3: movement towards one variable's center
    var_int_im = as.im(data.frame(x = simpoRes$quadsenv$X, y = simpoRes$quadsenv$Y, 
                                  z = simpoRes$quadsenv[, n_var]))
    
    var_rast = raster(var_int_im)
    
    vals.v1 <- values(var_rast)
    coord <- xyFromCell(var_rast, 1:ncell(var_rast))
    comb.v1 <- cbind(coord,vals.v1)
    
    C.v1 = as.data.frame(comb.v1)
    mc.x = C.v1$x[C.v1$vals.v1==max(C.v1$vals.v1)]
    mc.y = C.v1$y[C.v1$vals.v1==max(C.v1$vals.v1)]
    
    mc.var <- cbind(mc.x, mc.y)
    
    if(isTRUE(Allpt)){
      Sp1_y2_mov3 = simpoRes$PO
      Sp1_y2_mov3$x = rescale(Sp1_y2_mov3$x, to=c(0, 100-p_scale*(100-mc.x)))
      Sp1_y2_mov3$y = rescale(Sp1_y2_mov3$y, to=c(0, 100-p_scale*(100-mc.y)))

      par(mfrow=c(1,2))
      plot(simpoRes$PO, main="before", col="white")
      points(simpoRes$PO, pch=16, col=rgb(0.8,0,0, 0.4))
      
      plot(Sp1_y2_mov3, main="after", pch=16, cols=rgb(0,0.8,1, 0.4))
      points(mc.var, pch="*", col="red", cex=1.8)
      par(mfrow=c(1,1))
      
    }else{
      
      # choose which points to move
      Sp1_sub = sample(1:simpoRes$PO$n, floor(p_move*simpoRes$PO$n))
      Sp1_yr2_static  = simpoRes$PO[-Sp1_sub]
      
      Sp1_y2_mov3 = simpoRes$PO[Sp1_sub]
      Sp1_y2_mov3b = simpoRes$PO[Sp1_sub]
      Sp1_y2_mov3$x = rescale(Sp1_y2_mov3$x, to=c(0, 100-p_scale*(100-mc.x)))
      Sp1_y2_mov3$y = rescale(Sp1_y2_mov3$y, to=c(0, 100-p_scale*(100-mc.y)))
      
      Sp1_mov_yr2= superimpose(static=Sp1_yr2_static, move=Sp1_y2_mov3)
      
      par(mfrow=c(1,2))
      plot(simpoRes$PO, main="before", col="white")
      points(simpoRes$PO, pch=16, col=rgb(0.8,0,0, 0.4))
      points(Sp1_y2_mov3b, pch=16, col=rgb(0,0.8,1, 0.4))
      
      plot(Sp1_mov_yr2, main="after", pch=c(16,16), cols=c(rgb(0.8,0,0, 0.4), rgb(0,0.8,1, 0.4)), legend=F)
      points(mc.var, pch="*", col="red", cex=1.8)
      par(mfrow=c(1,1))
    }
    
    return(list(PO=simpoRes$PO, PO_mov=Sp1_y2_mov3, center = mc.var))
    
  }
}


# Functions from Renner et al., 2019
clogloginv = function(lin) {1 - exp(-exp(lin))}

newenv.var = function(sp.xy, env.grid, env.scale, coord = c("X","Y"), file.name = NA)
{
  convert = FALSE
  if (any(lapply(env.grid, class) == "factor"))
  {
    convert  = TRUE
    out.grid = CatConvert(env.grid)
    env.grid = out.grid$X
  }
  x.dat   = sp.xy[,which(names(sp.xy) == coord[1])]
  y.dat   = sp.xy[,which(names(sp.xy) == coord[2])]
  x.back  = env.grid[,which(names(env.grid) == coord[1])]
  y.back  = env.grid[,which(names(env.grid) == coord[2])]
  x.col   = which(names(env.grid) == coord[1])
  y.col   = which(names(env.grid) == coord[2])
  var.col = setdiff(1:dim(env.grid)[2], c(x.col, y.col))
  s.res   = SpatRes(env.grid)
  
  sp.dat        = as.data.frame(matrix(NA, length(x.dat), length(var.col)))
  names(sp.dat) = names(env.grid[var.col])
  
  for (var in 1:length(var.col))
  {
    loop.scale = min(c(s.res$x.step, s.res$y.step))
    loc        = which(is.na(sp.dat[,var]))
    while(sum(is.na(sp.dat[,var])) > 0)
    {
      loc = which(is.na(sp.dat[,var]))
      sp.dat[loc, var] = newinterp(sp.xy[loc,], loop.scale, env.grid[,var.col[var]], env.grid, coord = c("X","Y"))
      loop.scale = loop.scale*2
    }
    cat(paste("Calculating species environmental data for variable:", names(sp.dat)[var], "\n"))
    flush.console()
  }
  
  sp.dat = data.frame(x.dat, y.dat, sp.dat)
  names(sp.dat)[1:2] = c("X", "Y")
  if (is.na(file.name) == FALSE)
  {
    save.name = paste(file.name, ".RData", sep = "")
    save(sp.dat, file = save.name)
    print(paste("Output saved in the file", save.name))
  }
  if (convert == TRUE)
  {
    sp.dat = list(X = sp.dat, cat.names = out.grid$cat.names)
  }
  sp.dat
}

SpatRes = function(env.grid, coord = c("X", "Y"))
{
  x.col   = which(names(env.grid) == coord[1])
  y.col   = which(names(env.grid) == coord[2])
  x.uq    = sort(unique(env.grid[, x.col]))
  y.uq    = sort(unique(env.grid[, y.col]))
  n.dec   = max(unlist(lapply(x.uq, DecimalCount)))
  x.diff  = diff(x.uq)
  y.diff  = diff(y.uq)
  x.dec   = unlist(lapply(x.diff, DecimalCount))
  y.dec   = unlist(lapply(y.diff, DecimalCount))
  x.step  = min(floor(x.diff*10^max(x.dec) + 0.1))/(10^max(x.dec))
  y.step  = min(floor(y.diff*10^max(y.dec) + 0.1))/(10^max(y.dec))
  return(list(x.step = x.step, y.step = y.step))
}

DecimalCount = function(x, max.dec = max(10, max(nchar(x))), tol = 1.e-1)
{
  digits  = 0:max.dec
  x.diff = (x - round(x, digits))/(10^(-1*(digits + 3)))
  num.dec = digits[min(which(abs(x.diff) < tol))]
  num.dec
}

newinterp = function(sp.xy, sp.scale, f, back.xy, coord = c("X","Y"))
{
  options(scipen = 999)
  x.dat   = sp.xy[,which(names(sp.xy) == coord[1])]
  y.dat   = sp.xy[,which(names(sp.xy) == coord[2])]
  x.back  = back.xy[,which(names(back.xy) == coord[1])]
  y.back  = back.xy[,which(names(back.xy) == coord[2])]
  
  grid    = data.table(x.back, y.back, f, key = c("x.back", "y.back"))
  
  ux    = sort(unique(x.back))
  uy    = sort(unique(y.back))
  
  x.col   = which(names(back.xy) == coord[1])
  y.col   = which(names(back.xy) == coord[2])
  
  x.step  = ux[2] - ux[1]
  y.step  = uy[2] - uy[1]
  
  x.o = min(back.xy[,x.col]) - floor(min(back.xy[,x.col])/x.step)*x.step
  y.o = min(back.xy[,y.col]) - floor(min(back.xy[,y.col])/y.step)*y.step
  
  x.1   = floor((x.dat - x.o)/sp.scale)*sp.scale + x.o
  y.1   = floor((y.dat - y.o)/sp.scale)*sp.scale + y.o
  x.2   = pmin(x.1 + sp.scale, max(ux))
  y.2   = pmin(y.1 + sp.scale, max(uy))
  
  w11   = (x.2 - x.dat)*(y.2 - y.dat)/((x.2 - x.1)*(y.2 - y.1))
  w12   = (x.2 - x.dat)*(y.dat - y.1)/((x.2 - x.1)*(y.2 - y.1))
  w21   = (x.dat - x.1)*(y.2 - y.dat)/((x.2 - x.1)*(y.2 - y.1))
  w22   = (x.dat - x.1)*(y.dat - y.1)/((x.2 - x.1)*(y.2 - y.1))
  
  f11   = grid[list(x.1, y.1)]$f
  f12   = grid[list(x.1, y.2)]$f
  f21   = grid[list(x.2, y.1)]$f
  f22   = grid[list(x.2, y.2)]$f
  
  c11 = 1 - is.na(f11)
  c12 = 1 - is.na(f12)
  c21 = 1 - is.na(f21)
  c22 = 1 - is.na(f22)
  
  env.wt.mat = cbind(f11*w11*c11, f12*w12*c12, f21*w21*c21, f22*w22*c22) 
  
  f.interp = apply(env.wt.mat, 1, sum, na.rm = TRUE)/(w11*c11 + w12*c12 + w21*c21 + w22*c22)
  f.interp
}
