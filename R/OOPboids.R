# BOIDS inspired

# Physics based engine working with vectors
# Using basic physics engine tailored to R code that uses position, velocity, acceleration and force
# pos, vel, and acc are vectors that are defined as points x,y that represent vectors from the origin to x,y and hence have a magnitude and direction



# detect density of other nearby cells
# k-nearest density estimation (elife Ding et al. 2019)
# individual-level behavioural quantification

euclidean_dist <- function(x, y) sqrt(sum((x - y)^2))

# Functions needed to implement the boids algorithm
# each boid having a limited perception of its environment, results in emergent behavior of flocking

#steering force
# desired velocity - current velocity

#alignment,
#desiredVelocity is calculated by averaging the velocities of boids in radius
#steering itself is based on desired velocity minus current velocity (Craig Reynhold)
align <- function(obj, lsobj, perceptionRadius){
  steeringForce <- c(0,0)
  n <- 0 #total n other boids in perception radius
  for(i in 1:length(lsobj)){
    if(euclidean_dist(obj@pos, lsobj[[i]]@pos) < perceptionRadius){
      if(obj@id != lsobj[[i]]@id){ #non self
        steeringForce <- steeringForce + lsobj[[i]]@vel #adding all vectors
        n <- n + 1
      }
    }
  }
  if(n>0){
    steeringForce <- steeringForce / n
    steeringForce <- setMagnitude(steeringForce,obj@maxSpeed)
    steeringForce <- steeringForce - obj@vel
  }
  steeringForce <- setMagnitude(steeringForce,obj@maxForce)
  return(steeringForce)
}


#cohesion
cohesion <- function(obj, lsobj, perceptionRadius){
  steeringForce <- c(0,0)
  n <- 0 #total n other boids in perception radius
  for(i in 1:length(lsobj)){
    if(euclidean_dist(obj@pos, lsobj[[i]]@pos) < perceptionRadius){
      if(obj@id != lsobj[[i]]@id){ #non self
        steeringForce <- steeringForce + lsobj[[i]]@pos #adding all vectors
        n <- n + 1
      }
    }
  }
  if(n>0){
    steeringForce <- steeringForce / n
    steeringForce <- steeringForce - obj@pos
    steeringForce <- setMagnitude(steeringForce,obj@maxSpeed)
    steeringForce <- steeringForce - obj@vel
    steeringForce <- setMagnitude(steeringForce,obj@maxForce)
  }
  #steeringForce = setMagnitude(steeringForce,obj@maxForce)
  return(steeringForce)
}

#separation
separation <- function(obj, lsobj, perceptionRadius){
  steeringForce <- c(0,0)
  n <- 0 #total n other boids in perception radius
  for(i in 1:length(lsobj)){
    d <- euclidean_dist(obj@pos, lsobj[[i]]@pos)
    if(d < perceptionRadius){
      if(obj@id != lsobj[[i]]@id){ #non self
        #vector that points from the other to the obj's position
        #diff <- obj@pos - lsobj[[i]]@pos
        diff <-  obj@pos + lsobj[[i]]@pos
        #force needs to be inversely related to the distance between objects * (1/distance)
        diff <- diff / d
        steeringForce <- steeringForce + diff #adding all vectors
        n <- n + 1
      }
    }
  }
  if(n>0){
    steeringForce <- steeringForce / n
    steeringForce <- setMagnitude(steeringForce,obj@maxSpeed)
    steeringForce <- steeringForce - obj@vel
    steeringForce <- setMagnitude(steeringForce, obj@maxForce)
  }
  #steeringForce = setMagnitude(steeringForce,obj@maxForce)
  return(steeringForce)
}



#try with velocity vectors
# imagine vector from 0,0 to x,y representing a vector, with mangnitude (lenght) and diretion (angle to x=0 line
# make object with pos x y move at each timestep with a velocity v (with components v.x x.y )
# random vectors, ideally for velocity, perhaps we want a so called UNIT vector, that has a fixed magnitude of 1 and a random angle


unitVecor <- function(direction){
  alpha <- direction * pi / 180
  d <- 1  #Set d to 1 for unit vetors (magnitude=1)
  v.x <- d * cos(alpha)
  v.y <- d * sin(alpha)
  unitv <- c(v.x, v.y)
  unitv
}

#unitVecor(direction = runif(1,0,360))
#df <- c()
#for(i in 1:100){
#  df <- rbind(df, unitVecor(direction = runif(1,0,360)))
#}

#df %>% data.frame() %>% ggplot(aes(vx,vy))+geom_point()+coord_fixed()


#Calculate the magnitude of any given vector
#a2+b2=c2, c = sqr(a2+b2) (pythagoras)

#use magnitude to normalise vector to scale 1
#normalise a vector simply by dividing it by its magnitude (3,4) (magnitude =5) so normalised vector is (3/5, 4/5)
# so first calculate the magnitude and divide the vector by it
# set a magnitude of a vector to be fixed to particular none 1 value, one has to then multiply it
# perhaps also write a function to derive the direction of a given vector, this might be fun to record

#position, position of object relative to the origin 0,0
# velocity, #how position changes over time (each timeframe)
# acceleration, change of velocity over time (change direction or magnitude of the velocity vector)


#acceleration
# newtons low, Force=mass*acceleration

#function to calculate the magnitude of a given vector (as defined by x,y coordinates (relative to the origin))
magnitude <- function(vector){
  c <- sqrt((vector[1]^2+ vector[2]^2))
  c
}

#calculate the direction of a vector defined by x,y coordinates, relative to the x-axis (expressed in degrees)
# implementation uses atan2, the 2-argument arctangent. by which theta=atan2(y,x) is tha anlge (in radians with -pi < theta < pi)
# theta is expressed and returned as an angle measured in degrees (convered by theta_degree = theta_radians * 180/pi )
direction <-function(vector){
  theta <- atan2(vector[2],vector[1])
  theta * (180/pi)
}

#change a vector to have a set magnitude but retain the direction
# extracts the direction of a vector, and turns it into a unit vector of mangitude 1 before multiplying it by the desired magnitude
setMagnitude <- function(vector,magnitude){
  vector.o <- unitVecor(direction=direction(vector)) * magnitude
  vector.o
}


s4list2df <- function(lsobj){
   #faster implementation based on lists and data.frame function outside the for-loop
  tmpls <- list()
      for(i in names(attributes(lsobj[[1]]))[1:length(names(attributes(lsobj[[1]])))-1]){
      tmpls[[i]] <- unlist(lapply(lsobj, function(x) slot(x, i)))
      if(length(tmpls[[i]])>length(lsobj)){
        tmpls[[i]] <- cbind.data.frame(split(tmpls[[i]], rep(1:2, times=length(tmpls[[i]])/2)), stringsAsFactors=F)
      }
    }
  cnames <- c()
  for(q in  names(tmpls)){
    if(!is.null(dim(tmpls[[q]]))){
      p <- paste(q,1:dim(tmpls[[q]])[2],sep='')
    }else{
      p <- q
    }
    cnames <- c(cnames,p)
  }

  tmpdf <- data.frame(Reduce(cbind.data.frame, tmpls))
  colnames(tmpdf) <- cnames
  tmpdf
}

print.cell <- function(obj){
  for(n in names(obj)){
    cat(n,": " ,obj[[n]], "\n")
  }
}

setClass("cell",
         representation(id='character',
                        pos='numeric', #position vector (x and y)
                        vel='numeric', #velocity
                        acc='numeric', #acceleration
                        maxSpeed='numeric', #maximum speed
                        maxForce='numeric', #maximum force
                        r='numeric'), #radius of object
         prototype(id='cell1',
                   pos=c(x=0,y=0),
                   vel=c(x=0,y=1),
                   acc=c(0,0),
                   maxSpeed=4,
                   maxForce=0.01,
                   r=16))




# newtons 2nd law of motion
# Force = mass * acceleration
# assuming mass =1, then force=acceleration
# to implement Net force, or the sum of all forces acting on the object, we add force to the current acceleration (and reset acceleration each time)
# this allows a generic function to be called consecutively when different forces are applied (for example gravity and wind)
applyForce <- function(obj, force){
  obj@acc <- obj@acc + force
  obj
}


update <- function(obj){
  #limit acceleration
  if(magnitude(obj@acc)>1){
    obj@acc <- setMagnitude(obj@acc, magnitude =1)
  }

  #add acceleration vector to velocity vector
  obj@vel <- obj@vel + obj@acc

  #limit velocity
  if(magnitude(obj@vel) > obj@maxSpeed){
    obj@vel <- setMagnitude(obj@vel, magnitude = obj@maxSpeed)
  }

  #add velocity vector to possition
  obj@pos <- obj@pos + obj@vel #vector adding, position vector + velocity vector
  obj@acc <- c(0,0)
  obj
}


edges <- function(obj, limits=100){
  if(obj@pos[['x']] > limits){
    obj@pos[['x']]= -limits
  }else if(obj@pos[['x']] < -limits){
    obj@pos[['x']]=limits
  }
  if(obj@pos[['y']] > limits){
    obj@pos[['y']]= -limits
  }else if(obj@pos[['y']] < -limits){
    obj@pos[['y']]=limits
  }
  obj
}



v_limits <- 100
n_steps <- 1000



#generate multiple prey particles
generate_cells<- function(n_cell = 1,limits=100, positions='random'){
  ls.cell<-list()
  for(i in 1:n_cell){
    ls.cell[[i]] <- new("cell",
                        id=paste0('cell',i),
                        pos=c(x=runif(1,-limits,limits),y=runif(1,-limits,limits)),
                        vel=unitVecor(direction=runif(1,0,360)),
                        acc=c(0,0),
                        maxSpeed=1,
                        maxForce=0.005,
                        r=16)
  }
  ls.cell
}


#generate with random start positions
ls.cell <- generate_cells(n_cell = 30, limits=v_limits)

cell_tracker <- list()
for(i in 1:n_steps){
  snapshot.cell <- s4list2df(ls.cell)
  cell_tracker[[i]] <- snapshot.cell #rbind(cell_tracker, cbind('step'=i,snapshot.cell ))

  for(i_cell in 1:length(ls.cell)){
      #limits
      #ls.cell[[i_cell]] <- edges(ls.cell[[i_cell]], limits = v_limits)

      #accelrate towards origin - imaging a graviational pull
      #ls.cell[[i_cell]]@acc <- c(0,0) - ls.cell[[i_cell]]@pos

      # random walk?
      #force=c(runif(n=1,min=-2,max=2),runif(n=1,min=-2,max=2))
      #ls.cell[[i_cell]] <- applyForce(ls.cell[[i_cell]],force)

      #boids algorithm
      #align
      alignmentforce <- align(obj = ls.cell[[i_cell]], lsobj = ls.cell,perceptionRadius = 1 )
      ls.cell[[i_cell]] <- applyForce(ls.cell[[i_cell]], alignmentforce)

      #cohesion
      cohesionforce <- cohesion(obj = ls.cell[[i_cell]], lsobj = ls.cell,perceptionRadius = 100 )
      ls.cell[[i_cell]] <- applyForce(ls.cell[[i_cell]], cohesionforce)

      #separation
      #separationforce <- separation(obj = ls.cell[[i_cell]], lsobj = ls.cell,perceptionRadius = 5 )
      #ls.cell[[i_cell]] <- applyForce(ls.cell[[i_cell]], separationforce)

      #update
      ls.cell[[i_cell]] <- update(ls.cell[[i_cell]])
  }
}





cell_out <- dplyr::bind_rows(cell_tracker, .id = "step")
cell_out$step <- as.numeric(cell_out$step)


cell_out %>%
  ggplot2::ggplot(ggplot2::aes(pos1,pos2)) +
  ggplot2::geom_point() +
  gganimate::transition_time(step) +
  ggnewscale::new_scale_color() +
  ggplot2::scale_color_manual(values=c('grey','red','blue'))




cell_out %>%
  ggplot2::ggplot(ggplot2::aes(pos1,pos2)) +
  ggplot2::geom_path(aes(group=id,color=step),alpha=0.9) +
  gganimate::transition_reveal(along = step) +
  ggplot2::scale_color_gradientn(colours=c("#E4F00A", "white", "#22FF00"))+
  ggplot2::geom_point(aes(x=0,y=0),color='white')+
  theme(panel.background = element_rect(fill='black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())






cell_out %>%
  ggplot2::ggplot(ggplot2::aes(pos1,pos2)) +
  ggplot2::geom_path(aes(group=id,color=step),alpha=0.3) +
  gganimate::transition_reveal(along = step) +
  ggnewscale::new_scale_color() +
  ggplot2::scale_color_manual(values=c('grey','red','blue'))

