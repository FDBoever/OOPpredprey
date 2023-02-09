#predator prey in 2D space interception
#https://math.stackexchange.com/questions/143932/calculate-point-given-x-y-angle-and-distance
#Since π rad = 180°, degrees (d) can be converted to radians (r) using r = d · π ∕ 180 and the conversion of radians to degrees is d = r · 180 ∕ π.



# in one microliter, or one cubic mm
# or a grid of 1000x1000 - 1000 micrometer for 1 mm lenght



#define classes

a <- list(id='pred00001', #unique identifier
          size=5, #size as measured in cell diameter (micrometer)
          speed=1, #speed
          pos_x = 1, #position along x
          pos_y = 1 #position along y
          )

class(a) <- 'predator'


print.predator <- function(obj){
  for(n in names(obj)){
    cat(n,": " ,obj[[n]], "\n")
  }
}

print(a)


#define S4 class
setClass('predator',
         slots=list(id='character',
                    size='numeric',
                    speed='numeric',
                    pos_x='numeric',
                    pos_y='numeric'))


setClass("predator",
         representation(id='character',
                        size='numeric',
                        speed='numeric',
                        pos_x='numeric',
                        pos_y='numeric',
                        direction='numeric',
                        energy='numeric',
                        max_energy='numeric'), #in degrees
         prototype(id='pred1',
                   size=5,
                   speed=1,
                   pos_x=0,
                   pos_y=0,
                  direction=45,
                  energy=1000,
                  max_energy=1000))



setClass("prey",
         representation(id='character',
                        size='numeric',
                        speed='numeric',
                        pos_x='numeric',
                        pos_y='numeric',
                        direction='numeric',
                        status='character'), #in degrees
         prototype(id='pred1',
                   size=1,
                   speed=0.5,
                   pos_x=0,
                   pos_y=0,
                   direction=0,
                   status='alive'))



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

#setGeneric("pos_x", function(x) standardGeneric("pos_x"))
#setGeneric("pos_x<-", function(x, value) standardGeneric("pos_x<-"))

#setMethod("pos_x", "predator", function(x) x@pos_x)
#setMethod("pos_x<-", "predator", function(x, value) {
#  x@pos_x <- pos_x
#  x
#})


move <- function(obj, limits=NULL){
  # move point distance d over angle alpha (in radians)

  d <- obj@speed
  alpha <- obj@direction * pi / 180

  if(is.null(limits)){
    obj@pos_x = obj@pos_x + (d * cos(alpha))
    obj@pos_y = obj@pos_y + (d * sin(alpha))
    obj
    }
  else{
    while(obj@pos_x + (d * cos(alpha)) > limits | obj@pos_x + (d * cos(alpha)) < -limits  | obj@pos_y + (d * sin(alpha)) > limits  | obj@pos_y + (d * sin(alpha)) < -limits){
      alpha <- runif(1,0,360) * pi / 180
    }
    obj@pos_x = obj@pos_x + (d * cos(alpha))
    obj@pos_y = obj@pos_y + (d * sin(alpha))
    obj

  }
}


alter_energy <- function(obj, value){
  obj@energy = obj@energy + value
  if(obj@energy>obj@max_energy){
    obj@energy=obj@max_energy
  }
  obj
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



#library(dplyr)
#library(ggplot2)




#generate multiple prey particles
generate_prey<- function(n_prey = 100,limits=100){
  ls.prey<-list()
  for(i in 1:n_prey){
    ls.prey[[i]] <- new("prey", id=paste0('prey',i),size=1, speed=1,
                        pos_x=runif(1,-limits,limits),
                        pos_y=runif(1,-limits,limits),
                        direction=runif(1,0,360))
  }
  ls.prey
}



#extract dataframe from list of s4 objects
s4list2df <- function(lsobj){
  data.frame('pos_x'= unlist(lapply(lsobj, function(x) slot(x, 'pos_x'))),
             'pos_y'= unlist(lapply(lsobj, function(x) slot(x, 'pos_y'))),
  'status'= unlist(lapply(lsobj, function(x) slot(x, 'status'))),
              'size'= unlist(lapply(lsobj, function(x) slot(x, 'size'))))

}

tmp_prey_pos <- s4list2df(ls.prey)





#--------------------#
#--------------------#
#Create an object
s <- new("predator", id='pred00002',size=5, speed=1, pos_x=0,pos_y=0,direction=45, energy=1000)
show(s)


ls.prey <- generate_prey(n_prey = 1000, limits=250)

df <- c()
for(i in 1:1000){
  df <- rbind(df,
              c('step'=i,
                x=s@pos_x,
                y=s@pos_y,
                energy=s@energy))


  s <- move(s, limits=250)
  s <- alter_energy(s, value= -1)

  #tumble
  if(runif(1,0,1)>0.8){
    s@direction <- runif(1,0,360)
  }

  #check particle interception
  for(i_pr in 1:length(ls.prey)){
    if(ls.prey[[i_pr]]@status == 'alive'){
      if(circle.interception(x0=s@pos_x,
                             y0=s@pos_y,
                             R0=s@size/2,
                             x1=ls.prey[[i_pr]]@pos_x,
                             y1=ls.prey[[i_pr]]@pos_y,
                             R1=ls.prey[[i_pr]]@size/2
      )){
        message('Bingo!')

        #set prey to eaten
        ls.prey[[i_pr]]@status <- 'eaten'

        #increase_energy - full charge?
        s <- alter_energy(s,value=250)
      }
    }
  }

}

tmp_prey_pos <- s4list2df(ls.prey)

df <- data.frame(df)

df %>% ggplot(aes(x,y)) +
  geom_point(aes(color=energy)) +
  geom_path(aes(color=energy))+
  ggnewscale::new_scale_color() +
  geom_point(data=tmp_prey_pos,aes(pos_x,pos_y,color=status)) +
  scale_color_manual(values=c('blue','red'))




df %>% ggplot(aes(step,energy)) +
  geom_point(aes(color=energy)) +
  geom_path(aes(color=energy))+
  ggnewscale::new_scale_color() +
  geom_point(data=tmp_prey_pos,aes(pos_x,pos_y,color=status)) +
  scale_color_manual(values=c('blue','red'))





