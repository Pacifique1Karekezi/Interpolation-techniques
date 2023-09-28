# Load libraries
library(dplyr)
library(gstat)
library(sp)
library(ggplot2)
library(gridExtra)
library(INLA)
library(maptools)

# Funtion for suppressing output
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# Set seed
set.seed(10)

# Load data
data(meuse)
str(meuse)


# Show grid locations
data(meuse.grid)
ggplot(data = meuse.grid, aes(x, y)) +
  geom_point(size = 0.5, alpha = 0.5) + 
  labs(title = 'Meuse Grid') +
  coord_equal() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Create bubble plot
ggplot(data = meuse, aes(x, y)) +
  geom_point(aes(size = zinc), color = "blue", alpha = 3/4) + 
  labs(title = 'Zinc Concentration (ppm)') +
  coord_equal() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 8),
        legend.position = c(0.15, 0.85))

#############KRIGING#################

# Create SpatialPointsDataFrame and assign data its coordinate reference system
coordinates(meuse) <- ~x+y
proj4string(meuse) <- CRS("+init=epsg:28992")

# Perform similar operation on the grid
coordinates(meuse.grid) = ~x+y
proj4string(meuse.grid) <- CRS("+init=epsg:28992")
gridded(meuse.grid) = TRUE

# Fit variogram
vgm <- variogram(log(zinc) ~ dist, meuse)
fit.vgm <- fit.variogram(vgm, vgm("Sph"))

# Universal kriging
krg <- quiet(krige(log(zinc) ~ dist, meuse, meuse.grid, model = fit.vgm))

# Add estimates to meuse.grid
meuse.grid$zinc.krg <- krg$var1.pred
meuse.grid$zinc.krg.sd <- sqrt(krg$var1.var)

# Create plot of kriging results
meuse_spdf <- as.data.frame(meuse.grid)
ggplot(meuse_spdf, aes(x = x, y = y, fill = zinc.krg)) +
  geom_tile() +
  scale_fill_distiller('Zinc', palette = 'Spectral',
                       na.value = 'transparent',
                       limits = c(4, 7.2)) +
  labs(title = 'Log-Concentration of Zinc (ppm)',
       subtitle = 'Universal Kriging') +
  coord_equal() +
  theme(plot.title = element_text(margin = margin(b = 2), size = 12,
                                  hjust = 0.0, color = 'black',
                                  face = quote(bold)),
        plot.subtitle = element_text(margin = margin(b = 2), size = 10,
                                     hjust = 0.0, color = 'black',
                                     face = quote(italic)),
        line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank())

######INLA##########

# Boundary
meuse.bdy <- unionSpatialPolygons(
  as(meuse.grid, "SpatialPolygons"), rep (1, length(meuse.grid))
)

# Define mesh
coo <- meuse@coords
pts <- meuse.bdy@polygons[[1]]@Polygons[[1]]@coords
mesh <- inla.mesh.2d(loc.domain = pts, max.edge = c(150, 500),
                     offset = c(100, 250))

par(mar = c(0, 0, 0.8, 0))
plot(mesh, asp = 1)
points(coo, col = 'red')
lines(pts, col = 'blue', lwd = 3)

# Create SPDE
meuse.spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
A.meuse <- inla.spde.make.A(mesh = mesh, loc = coordinates(meuse))
s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = meuse.spde$n.spde)

# Create data structure
meuse.stack <- inla.stack(data  = list(zinc = meuse$zinc),
                          A = list(A.meuse, 1),
                          effects = list(c(s.index, list(Intercept = 1)),
                                         list(dist = meuse$dist)),
                          tag = "meuse.data")

# Create data structure for prediction
A.pred <- inla.spde.make.A(mesh = mesh, loc = coordinates(meuse.grid))
meuse.stack.pred <- inla.stack(data = list(zinc = NA),
                               A = list(A.pred, 1),
                               effects = list(c(s.index, list (Intercept = 1)),
                                              list(dist = meuse.grid$dist)),
                               tag = "meuse.pred")

# Join stack
join.stack <- inla.stack(meuse.stack, meuse.stack.pred)

# Fit model
form <- log(zinc) ~ -1 + Intercept + dist + f(spatial.field, model = spde)

m1 <- inla(form, data = inla.stack.data(join.stack, spde = meuse.spde),
           family = "gaussian",
           control.predictor = list(A = inla.stack.A(join.stack), compute = TRUE),
           control.compute = list(cpo = TRUE, dic = TRUE))

# Summary of results
summary(m1)

# Get predicted data on grid
index.pred <- inla.stack.index(join.stack, "meuse.pred")$data
meuse.grid$zinc.spde <- m1$summary.fitted.values[index.pred, "mean"]
meuse.grid$zinc.spde.sd <- m1$summary.fitted.values[index.pred, "sd"]

# Plot posterior means using ggplot
meuse_spdf <- as.data.frame(meuse.grid)
p1 <- ggplot(meuse_spdf, aes(x = x, y = y, fill = zinc.krg)) +
  geom_tile() +
  scale_fill_distiller('Zinc', palette = 'Spectral',
                       na.value = 'transparent',
                       limits = c(4, 7.2)) +
  labs(title = 'Log-Concentration of Zinc (ppm)',
       subtitle = 'Universal Kriging') +
  coord_equal() +
  theme(plot.title = element_text(margin = margin(b = 2), size = 12,
                                  hjust = 0.0, color = 'black',
                                  face = quote(bold)),
        plot.subtitle = element_text(margin = margin(b = 2), size = 10,
                                     hjust = 0.0, color = 'black',
                                     face = quote(italic)),
        line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank())

p2 <- ggplot(meuse_spdf, aes(x = x, y = y, fill = zinc.spde)) +
  geom_tile() +
  scale_fill_distiller('Zinc', palette = 'Spectral',
                       na.value = 'transparent',
                       limits = c(4, 7.2)) +
  labs(title = 'Log-Concentration of Zinc (ppm)',
       subtitle = 'INLA-SPDE') +
  coord_equal() +
  theme(plot.title = element_text(margin = margin(b = 2), size = 12,
                                  hjust = 0.0, color = 'black',
                                  face = quote(bold)),
        plot.subtitle = element_text(margin = margin(b = 2), size = 10,
                                     hjust = 0.0, color = 'black',
                                     face = quote(italic)),
        line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank())

grid.arrange(p1, p2, nrow = 1)


################################################
######KIGALI INLA################################
#####################################################
Kigalli<- st_read("kigali.shp")

kigali.bdy <- unionSpatialPolygons(
  as(as(Kigali,"Spatial"), "SpatialPolygons"), rep (1, length(as(Kigali,"Spatial")))
)


# Define mesh
coo <-simulated@coords
pts <- kigali.bdy@polygons[[1]]@Polygons[[1]]@coords
mesh <- inla.mesh.2d(loc.domain = pts, max.edge = c(150, 500),
                     offset = c(100, 250))

#plot(mesh)

par(mar = c(0, 0, 0.8, 0))
plot(mesh, asp = 1)
points(coo, col = 'red')
lines(pts, col = 'blue', lwd = 3)

# Create SPDE
simulated.spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
A.simulated <- inla.spde.make.A(mesh = mesh, loc = coordinates(simulated))
s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = simulated.spde$n.spde)

# Create data structure for prediction
A.pred <- inla.spde.make.A(mesh = mesh, loc = coordinates(as(Kigali,"Spatial")))
meuse.stack.pred <- inla.stack(data = list(sim1 = NA),
                               A = list(A.pred, 1),
                               effects = list(c(s.index, list (Intercept = 1))),
                               tag = "meuse.pred")


options(repos = c(getOption("repos"),
                  INLA="https://inla.r-inla-download.org/R/stable"))

# Install INLA and dependencies (from CRAN)
install.packages("INLA", dep = TRUE)



########################################################
##########INLA model##############################################
############################################################


##Simulated data in R

N <- 100 # 500, 5000, 25000, 100000
x <- rnorm(N, mean = 6, sd = 2)
y <- rnorm(N, mean = x, sd = 1)
data <- list(x = x, y = y, N = N)

##JAGS code
library(JAGS)

model <- function() {
  for(i in 1:N) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta * x[i]
  }
  alpha ~ dnorm(0, 0.001)
  beta  ~ dnorm(0, 0.001)
  tau   ~ dgamma(0.01, 0.01)
}
params <- c("alpha", "beta", "tau", "mu")
jags(
  data = data,
  param = params,
  n.chains = 3,
  n.iter = 50000,
  n.burnin = 5000,
  model.file = model
)

##INLA codes

res<-inla(y ~ x,
     family = "gaussian",
     data = data,
     control.predictor = list(link = 1)
)

summary(res)
