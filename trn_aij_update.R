updated_tg_aij <- function(tf, tg, H, L, w_m, w_trn) {
  tg_updated = 0
  if (tf <= L) {
    tg_updated = tg
  } else if (tf >= H) {
    val = 1/(w_m+w_trn)*(w_m*tg+w_trn*tf)
    if (val > tg) {
      tg_updated = val
    } else {
      tg_updated = tg
    }
  } else {
    low = updated_tg_aij(L,tg,H,L,w_m,w_trn)
    high = updated_tg_aij(H,tg,H,L,w_m,w_trn)
    #val = tg+(tf-L)*(high-low)
    m = (high-low)/(H-L)
    val = (tf-L)*m+tg
    if (val > tg) {
      tg_updated = val
    } else {
      tg_updated = tg
    }
  }
  return(tg_updated)
}

updated_tg_aij_vector <- function(tf, tg, H, L, w_m, w_trn) {
  return(lapply(tf, FUN=function(x) updated_tg_aij(x,tg,H,L,w_m,w_trn)))
}






H = 0.8
L = 0.2
w_m = 0.5
w_trn = 1
lim = c(0,1)
colors = c("red","orange","yellow","green","blue","purple","darkred","darkorange","darkgoldenrod","darkgreen","darkblue")
par(mfrow=c(1,1))
plot(c(),c(), xlim=lim, ylim=lim, ylab="a*_gj", xlab="a_fj", main="New TG aijs vs. TF aijs for Various TG aijs")
for (tg in seq(0,1,0.1)) {
  curve(updated_tg_aij_vector(x,tg,H,L,w_m,w_trn), add=TRUE, col=colors[as.integer(10*tg+1)])
  legend("bottomright", c("a_gj=0","a_gj=0.1","a_gj=0.2","a_gj=0.3","a_gj=0.4","a_gj=0.5","a_gj=0.6","a_gj=0.7","a_gj=0.8","a_gj=0.9","a_gj=1"), col=colors, pch=c(20), pt.cex=1,cex=0.8)
} 
