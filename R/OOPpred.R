
library(dplyr)
library(ggplot2)


#https://sites.google.com/site/shsgalbraith/non-class-pages/science-simulations/predator-prey-simulation


#predator prey in 2D space interception
#https://math.stackexchange.com/questions/143932/calculate-point-given-x-y-angle-and-distance
#Since π rad = 180°, degrees (d) can be converted to radians (r) using r = d · π ∕ 180 and the conversion of radians to degrees is d = r · 180 ∕ π.


#implementation in 3D space
#https://stackoverflow.com/questions/52781607/3d-point-from-two-angles-and-a-distance

#inspiring code adventures
#https://www.youtube.com/watch?v=X-iSQQgOd1A


#-----
# Currently working in confined space, boundaries

#------------------------------------------------#
# metabolic cost - derived from their size and seed
# metabolic cost 'energy used to stay alive/move per time unit' <- speed/2*size
# traits can be stored in tables, which can be modified by mutation, small number changes in trait values, that can affect behavior of cells and offspring

#look into NeuroEvolution of Augmenting Topologies (NEAT evolutionary algorithm)
# Stanley (2002), Evolving Neural Networks through augmenting topologies MIT Press journals
# RT-neat

# in one microliter, or one cubic mm
# or a grid of 1000x1000 - 1000 micrometer for 1 mm lenght



#define classes

print.predator <- function(obj){
  for(n in names(obj)){
    cat(n,": " ,obj[[n]], "\n")
  }
}


#define S4 class
setClass("predator",
         representation(id='character',
                        size='numeric',
                        speed='numeric',
                        pos_x='numeric',
                        pos_y='numeric',
                        direction='numeric',
                        energy='numeric',
                        max_energy='numeric',
                        generation='numeric',
                        status='character'), #in degrees
         prototype(id='pred1',
                   size=5,
                   speed=1,
                   pos_x=0,
                   pos_y=0,
                  direction=45,
                  energy=1000,
                  max_energy=1000,
                  generation=1,
                  status='alive'))

#generate multiple prey particles
generate_predators<- function(n_pred = 100,limits=100){
  ls.pred<-list()
  for(i in 1:n_pred){
    ls.pred[[i]] <- new("predator",
                        id=paste0('predator',i),
                        size=1,
                        speed=1,
                        pos_x=runif(1,-limits,limits),
                        pos_y=runif(1,-limits,limits),
                        direction=runif(1,0,360),
                        energy=1000,
                        max_energy=1000,
                        generation=1,
                        status='alive')
  }
  ls.pred
}


setClass("prey",
         representation(id='character',
                        size='numeric',
                        speed='numeric',
                        pos_x='numeric',
                        pos_y='numeric',
                        direction='numeric',
                        status='character',
                        motile='logical',
                        generation='numeric',
                        energy='numeric',
                        max_energy='numeric'),
         prototype(id='pred1',
                   size=1,
                   speed=0.25,
                   pos_x=0,
                   pos_y=0,
                   direction=0,
                   status='alive',
                   motile=TRUE,
                   generation=1,
                   energy=1000,
                   max_energy=1000))


#generate multiple prey particles
generate_prey<- function(n_prey = 100,limits=100){
  ls.prey<-list()
  for(i in 1:n_prey){
    ls.prey[[i]] <- new("prey",
                        id=paste0('prey',i),
                        size=1,
                        speed=0.25,
                        pos_x=runif(1,-limits,limits),
                        pos_y=runif(1,-limits,limits),
                        direction=runif(1,0,360),
                        status='alive',
                        motile=TRUE,
                        generation=1,
                        energy=1000,
                        max_energy=1000)
  }
  ls.prey
}

# Defining a function to display object details (setmethod will overwrite existing fcuntion for the specified object)
setMethod('show', 'predator',
          function(object){
            for(n in names(attributes(object))){
              cat(n,": " ,attr(object,n), "\n")
            }
          })

setMethod('show', 'prey',
          function(object){
            for(n in names(attributes(object))){
              cat(n,": " ,attr(object,n), "\n")
            }
          })


#leave limits == NULL when movement is unbound
# set litims (-limit,+limit) to move in confined space
move <- function(obj, limits=NULL, tumble_freq=0){

  # move point distance d over angle alpha (in radians)
  d <- obj@speed
  alpha <- obj@direction * pi / 180

  if(is.null(limits)){
    obj@pos_x = obj@pos_x + (d * cos(alpha))
    obj@pos_y = obj@pos_y + (d * sin(alpha))
    }
  else{
    while(obj@pos_x + (d * cos(alpha)) > limits | obj@pos_x + (d * cos(alpha)) < -limits  | obj@pos_y + (d * sin(alpha)) > limits  | obj@pos_y + (d * sin(alpha)) < -limits){
      alpha <- runif(1,0,360) * pi / 180
    }
    obj@pos_x = obj@pos_x + (d * cos(alpha))
    obj@pos_y = obj@pos_y + (d * sin(alpha))
  }

  #tumble
  if(tumble_freq>0){
    if(runif(1,0,1)<=tumble_freq){
      obj@direction <- runif(1,0,360)
    }
  }
  obj
}


alter_energy <- function(obj, value){
  obj@energy = obj@energy + value
  if(obj@energy>obj@max_energy){
    obj@energy=obj@max_energy
  }
  obj
}



divide_org <- function(lsobj, i, K = NULL){
  if(!is.null(K)){
    if(K>length(lsobj[unlist(lapply(lsobj, function(x) slot(x, 'status'))) == 'alive'])){
      lsobj <- append(lsobj, lsobj[[i]])
      lsobj <- append(lsobj, lsobj[[i]])
      lsobj[[length(lsobj)-1]]@generation <- lsobj[[i]]@generation + 1
      lsobj[[length(lsobj)]]@generation <- lsobj[[i]]@generation + 1
      lsobj[[length(lsobj)-1]]@energy <- lsobj[[i]]@max_energy
      lsobj[[length(lsobj)]]@energy <- lsobj[[i]]@max_energy
      #change direction of second offspring cell
      lsobj[[length(lsobj)]]@direction <- 360 - lsobj[[i]]@direction
      lsobj[[i]]@status <- 'divided'
    }
  }else{
    lsobj <- append(lsobj, lsobj[[i]])
    lsobj <- append(lsobj, lsobj[[i]])
    lsobj[[length(lsobj)-1]]@generation <- lsobj[[i]]@generation + 1
    lsobj[[length(lsobj)]]@generation <- lsobj[[i]]@generation + 1
    lsobj[[length(lsobj)-1]]@energy <- lsobj[[i]]@max_energy
    lsobj[[length(lsobj)]]@energy <- lsobj[[i]]@max_energy
    #change direction of second offspring cell
    lsobj[[length(lsobj)]]@direction <- 360 - lsobj[[i]]@direction
    lsobj[[i]]@status <- 'divided'
  }


  lsobj
}






#checking point interception
#checking if point x,y is in the circle with radius and centre
#returns true upon particle interception, target particle x,y
# center x,y and radius reflect position of predator with given size
point.interception <- function(x_target,y_target,x_center,y_center,radius){
  ifelse((x_target - x_center)^2 + (y_target - y_center)^2 < radius^2, TRUE, FALSE)
}

point.interception(10,0,0,0,5)


#Two circles intersect if, and only if, the distance between their centers is between the sum and the difference of their radii. Given two circles (x0, y0, R0) and (x1, y1, R1), the formula is as follows:
circle.interception <- function(x0,y0,R0, x1,y1,R1){
  ifelse((R0 - R1)^2 <= (x0 - x1)^2 + (y0 - y1)^2 & (x0 - x1)^2 + (y0 - y1)^2 <= (R0 + R1)^2,
         TRUE,
         FALSE)
}

circle.interception(0,0,5,0,0,5)










#extract dataframe from list of s4 objects
s4list2df <- function(lsobj){
  #tmpdf <- c()
  #for(i in names(attributes(lsobj[[1]]))[1:length(names(attributes(lsobj[[1]])))-1]){
  #  tmpdf <- data.frame(cbind(tmpdf, unlist(lapply(lsobj, function(x) slot(x, i)))))
  #}
  #tmpdf <- data.frame(tmpdf)
  #colnames(tmpdf) <- names(attributes(lsobj[[1]]))[1:length(names(attributes(lsobj[[1]])))-1]

  #faster implementation based on lists and data.frame function outside the for-loop
  tmpls <- list()
  for(i in names(attributes(lsobj[[1]]))[1:length(names(attributes(lsobj[[1]])))-1]){
    tmpls[[i]] <- unlist(lapply(lsobj, function(x) slot(x, i)))
  }
  tmpdf <- data.frame(Reduce(cbind.data.frame, tmpls))
  colnames(tmpdf) <- names(tmpls)
  tmpdf
}


system.time({
  tmp_prey_pos <- s4list2df(ls.prey)
})



#--------------------#
#--------------------#
#Create an object
s <- new("predator", id='pred00002',size=5, speed=1, pos_x=0,pos_y=0,direction=45, energy=1000)
show(s)


#each step of this iteration represent an arbiraty lenght of time
#the exact lenght of these steps in terms of real time is not defined, although some aspects of the cells' lives are configured in terms of it
# for instance the speed by which they move (or (howlong they live, or how old they have to be in order to be able to give birth))


# prey population increases rapidly, while predator population may initially decrease slightly - lets see if we can reporduce Lotka-Volterra like oscillations - where the peak in predator numbers lags slighly behind the peak in fox numbers
# For lotkavolterra like dynamics we can consider using following attributes
# breeding-age=5, max-age=40, breeding probability=0.12, max_litter_size=2 (maximum number of cells they give birth to),
# age of cell, boolean whether it is alive, and location in the field
# move - may need to check whether there is space adjacently to where it is, easier in a grid-based space
# multiply if old enough and propperly fed
# eihter prey particles or predators should have a value that includes the food_value (food value associated with eating a prey particle - could be expressed in terms of energy) - energy in fact relates to the number of steps the predator can survive between prey interceptions
# predator enery keeps track of whether the predator is on the brink of starvation or not

#preys energy might want to start at 0 and increase over time (with the idea of taking up nutrients), when reaching threshold they could divide, and offspring set to 0 again

#need to implement assexual reproduction by division of one cell into 2, and settin energy levels



ls.prey <- generate_prey(n_prey = 30, limits=250)
ls.pred <- generate_predators(n_pred = 3, limits=250)


#set boundary
v_limits <- 250
v_limits <- 125
v_limits <- 50

K_prey <- 100
K_predator <- 10

n_steps <- 5000

ls.prey <- generate_prey(n_prey = 10, limits=v_limits)
ls.pred <- generate_predators(n_pred = 1, limits=v_limits)

#profvis::profvis({
predator_tracker <- list()
prey_tracker <- list()
for(i in 1:n_steps){
  snapshot.pred <- s4list2df(ls.pred)
  snapshot.prey <- s4list2df(ls.prey)
  predator_tracker[[i]] <- snapshot.pred #rbind(predator_tracker, cbind('step'=i,snapshot.pred ))
  prey_tracker[[i]] <- snapshot.prey #rbind(prey_tracker, cbind('step'=i, snapshot.prey))

  for(i_pred in 1:length(ls.pred)){
    if(ls.pred[[i_pred]]@status == 'alive'){

      #move predator
      ls.pred[[i_pred]] <- move(ls.pred[[i_pred]], limits=v_limits, tumble_freq=0.2)

      #cost of moving
      ls.pred[[i_pred]] <- alter_energy(ls.pred[[i_pred]], value= -1)

      #kill predator if running out of energy
      if(ls.pred[[i_pred]]@energy == 0){
        ls.pred[[i_pred]]@status ='dead'
      }


      if (i%%1000 == 0){
        ls.pred <- divide_org(lsobj = ls.pred, i = i_pred, K= K_predator)
      }

      #check prey interception
      # use snapshop to filter close proximity to predator
      neighbours <- as.numeric(
        rownames(
          snapshot.prey[snapshot.prey$pos_x > ls.pred[[i_pred]]@pos_x - (ls.pred[[i_pred]]@size*2) &
                          snapshot.prey$pos_x < ls.pred[[i_pred]]@pos_x + (ls.pred[[i_pred]]@size*2) &
                          snapshot.prey$pos_y > ls.pred[[i_pred]]@pos_y - (ls.pred[[i_pred]]@size*2) &
                          snapshot.prey$pos_y < ls.pred[[i_pred]]@pos_y + (ls.pred[[i_pred]]@size*2),]))
      if(length(neighbours)>0){
        for(i_pr in 1:length(ls.prey)){
          if(ls.prey[[i_pr]]@status == 'alive'){
            if(circle.interception(x0=ls.pred[[i_pred]]@pos_x,
                                   y0=ls.pred[[i_pred]]@pos_y,
                                   R0=ls.pred[[i_pred]]@size/2,
                                   x1=ls.prey[[i_pr]]@pos_x,
                                   y1=ls.prey[[i_pr]]@pos_y,
                                   R1=ls.prey[[i_pr]]@size/2
            )){
              message('prey interception')

              #set prey to eaten
              ls.prey[[i_pr]]@status <- 'eaten'

              #increase_energy - full charge?
              ls.pred[[i_pred]] <- alter_energy(ls.pred[[i_pred]],value=500)
            }
          }
      }
      }



    }
  }


  #---prey
  for(i_pr in 1:length(ls.prey)){
    if(ls.prey[[i_pr]]@status == 'alive'){
      #move prey
      if(ls.prey[[i_pr]]@motile == TRUE){
        ls.prey[[i_pr]] <- move(ls.prey[[i_pr]], limits=v_limits, tumble_freq = 0.2)
      }

      #divide at fixed time
      if (i%%250 == 0){
        ls.prey <- divide_org(lsobj = ls.prey, i = i_pr, K= K_prey)
      }
    }




  }
}

predator_out <- dplyr::bind_rows(predator_tracker, .id = "step")
predator_out$step <- as.numeric(predator_out$step)

prey_out <- dplyr::bind_rows(prey_tracker, .id = "step")
prey_out$step <- as.numeric(prey_out$step)
#})





#tmp_prey_pos <- s4list2df(ls.prey)
#predator_tracker <- data.frame(predator_tracker)
#predator_tracker$step <- as.integer(predator_tracker$step)


#p1 <- predator_tracker %>% ggplot(aes(x,y)) +
#  geom_point(aes(color=energy)) +
#  geom_path(aes(color=energy))+
#  ggnewscale::new_scale_color() +
#  geom_point(data=tmp_prey_pos,aes(pos_x,pos_y,color=status),size=0.5) +
#  scale_color_manual(values=c('blue','red'))+
#  theme(aspect.ratio=1)

#ggplot2::ggsave(p1,filename = 'bactoverlay_ggplot_static2.pdf',width=5,height = 5 )



#prey_tracker %>% ggplot(aes(pos_x,pos_y)) +
#  geom_point(aes(color=status)) +
#  ggnewscale::new_scale_color() +
#  gganimate::transition_time(step) +
#    theme(aspect.ratio=1)



pg <- predator_out %>%
  dplyr::filter(status=='alive') %>%
  ggplot2::ggplot(aes(pos_x,pos_y)) +
  ggplot2::geom_point(aes(color=energy)) +
  gganimate::transition_time(step) + ggnewscale::new_scale_color() +
  geom_point(data=prey_out %>% filter(status!='divided'),aes(pos_x,pos_y,color=status),size=0.5) +
  scale_color_manual(values=c('grey','red','blue'))

gganimate::animate(pg,#
                   fps=10)
                   #duration = 20)


datetime <- Sys.time()
gganimate::anim_save(filename=paste0('bactoverlay_gganimate_klo2o',datetime,'.gif'),
                     animation=gganimate::last_animation())







predator_out %>% filter(status=='alive') %>%
  group_by(step) %>%
  tally() %>%
  ggplot(aes(x=step, y=n))+geom_point()

prey_out %>% filter(status=='alive') %>%
  group_by(step) %>%
  tally() %>% ggplot(aes(x=step, y=n))+geom_point()




predator_tracker %>% ggplot(aes(x,y)) +
  geom_point(aes(color=energy)) +
  #geom_path(aes(color=energy))+
  gganimate::transition_time(step)


gganimate::anim_save(filename='test_gganimate.gif',animation=gganimate::last_animation())





predator_tracker %>% ggplot(aes(step,energy)) +
  geom_point(aes(color=energy)) +
  geom_path(aes(color=energy))+
  ggnewscale::new_scale_color() +
  geom_point(data=tmp_prey_pos,aes(pos_x,pos_y,color=status)) +
  scale_color_manual(values=c('blue','red'))


#predation - interesting form of interaction between organisms
#bear killing whatever although they gain a lot of their nutrition from bug, nuts, betties, honey
# evoltuioneay presures of hunting and being hunted#
# speed of prey, claws and teeth
#effect of predator-prey interaction that results in the mind-boggling amount of diversity we see in all ecosystems, from savanna to rockpools

#kill for energy
# flow of energy through nature, need energy to survive and reproduce
# grazing, organism eats other organism
# prey adaptations could be detection, capture and handling


