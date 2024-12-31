# Create samples of the recharge sites 
# Simulate 6 sites in the same AMIGO scenario

library(lhs)
library(data.table)
library(ggplot2)

# Defined constants
scenarios = 500/.7

## Constants from the IDF files
Xmin = 204475
Xmax = 241275

Ymin = 439250
Ymax = 460050

Amin = 100*100 # 1 hectare, decided with Perry and Ruud (07/02/2022)
Amax = 1000*1000 # 1sq km, decided with Perry and Ruud (07/02/2022)

Rmin = 5
Rmax = 25

spreading = 750
min_dist = 2*6500 # Extent>1cm from the steady state runs

## print the extent of the catchment away from the constant head boundary
dist = spreading*3
paste(Xmin+dist, Ymin+dist, Xmax-dist, Ymax-dist, sep = ",")

# Orthogonal array based Latin hypercube sample for the center co-ordinates
## Adds latin hypercube benifits to orthogonally distributed samples
## Only center co-ordinates as for the other (recharge and Area) 
## equal representation is very important
## Read para starting with "This problem was addressed by Tang (1993)" in 
## https://www.sciencedirect.com/science/article/pii/S0304407603002124

## Find the number of samples that can fit
Xrange = Xmax-Xmin # 3 x spreading length from left and right
Yrange = Ymax-Ymin

# We can simulate multiple at the same time. 
# Spliting the orthogonal grid into equal sections to be simulated together
X_splits = round(Xrange/min_dist)
Y_splits = round(Yrange/min_dist)
scenario_split = scenarios/(X_splits*Y_splits)
range_x_splits = Xrange/X_splits
range_y_splits = Yrange/Y_splits

# Nx = round(sqrt(scenario_split*range_x_splits/range_y_splits))
# Ny = round(scenario_split/Nx)

Nx = round(sqrt(scenarios*Xrange/Yrange))
Ny = round(scenarios/Nx)


## Create the orthogonal array
side = sqrt(Amax)
orthoGrid = expand.grid(x = seq(Xmin+side/2, Xmax-side/2, length.out = Nx),
                        y = seq(Ymin+side/2, Ymax-side/2, length.out = Ny))
setDT(orthoGrid)

## Scale the x and y values
setnames(orthoGrid, c("x", "y"))

orthoGrid[,rch_site:=0]
groups = 1
for (x_min_limit in seq(Xmin, Xmax, length.out = X_splits+1)[-1]){
  for (y_min_limit in seq(Ymin, Ymax, length.out = Y_splits+1)[-1]){
    
    orthoGrid[x<x_min_limit & y<y_min_limit & rch_site==0, rch_site:=groups]
    groups = groups+1
  }
}

orthoGrid[, scenario := 1:.N, rch_site]

Grid2 = randomLHS(orthoGrid[,.N], 2)

## Scaling
Grid2 = as.data.table(Grid2)
setnames(Grid2, c("recharge", "area"))

Grid2[, ":="(recharge = Rmin + recharge*(Rmax-Rmin),# Mean WTP for NL ~15500m3/day
             area = Amin + area*(Amax-Amin))]

orthoGrid = cbind(orthoGrid, Grid2)

orthoGrid[, ":="(side = sqrt(area))]
orthoGrid[, ":="(W_side = x - side/2,
                 E_side = x + side/2,
                 N_side = y + side/2,
                 S_side = y - side/2)]
library(dplyr)
windows(Xrange, Yrange)
orthoGrid%>%
  # filter(scenario==1)%>%
  ggplot()+
  geom_rect(aes(xmin=W_side, xmax=E_side, ymin=S_side, ymax=N_side,fill = as.factor(rch_site)))+
  coord_cartesian(xlim=c(Xmin, Xmax), ylim = c(Ymin, Ymax))

distances = orthoGrid[, .(distance = min(dist(.SD[,.(x, y)]))), .(scenario)]


fwrite(orthoGrid, "ortho_grid.csv")
