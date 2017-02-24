### R code from vignette source 'article.Rnw'

###################################################
### code chunk number 1: convolution-snippet
###################################################
K <- 1
## q <- 2*K
q <- 2*K+1
Kp1 <- K+1
even.deg <- q%%2 == 0

for(j in 1:q){
  print(bquote(g1[.(as.double(j))] <-
               inprod(bn[.(max(0,j-Kp1)+1):.(min(Kp1,j)),1],
                      bnrev[.(max(0,Kp1-j)+1):.(min(q+even.deg-j,K)+1),1])))
  print(bquote(g2[.(as.double(j))] <-
               inprod(bn[.(max(0,j-Kp1)+1):.(min(Kp1,j)),2],
                      bnrev[.(max(0,Kp1-j)+1):.(min(q+even.deg-j,K)+1),2])))
}


