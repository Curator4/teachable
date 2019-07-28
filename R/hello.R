# Made by Curator
# dev ready

C <- read.csv('caData.csv', header = FALSE)
labelx <- read.csv('LabelsCremas.csv', header = FALSE)
Zsupsred <- read.csv('Zsupsred.csv', header = FALSE)
Lsupsred <- read.csv('Lsupsred.csv', header = FALSE)

CAFc <- function(F, k, SP, LabelsSup, RorC) {

  minimum = min(min(F))

  if (minimum != 0) {
    F = F - minimum
  }

  # Positive attribute values
  F1 <- F

  # Negative attribute values
  F2 <- max(max(F1)) - F1

  # Doubled attribute matrix
  Fc <- cbind(F1, F2)


  s <- sum(Fc)
  cat <- dim(Fc)

  # The anylysis ----

  P <- Fc * (1 / s)

  svd <- svd(P)


  r <- as.matrix(P) %*% matrix(1, nrow = cat[2], ncol = 1)
  c <- as.matrix(t(P)) %*% matrix(1, nrow = cat[1], ncol = 1)

  Dr <- diag(as.vector(r))
  Dc <- diag(as.vector(c))
  Drh <- Dr ^0.5
  Dch <- Dc ^0.5

  svd <- svd(solve(Drh) %*% (as.matrix(P) - r %*% t(c)) %*% solve(Dch), nu = 380, nv = 24)

  rang <- qr(as.matrix(P) - r %*% t(c))$rank



  l <- diag(svd$d)[1:rang, 1:rang]
  U <- svd$u[, 1:rang]
  V <- svd$v[, 1:rang]

  Y <- solve(Drh) %*% U
  X <- solve(Dch) %*% V
  Yk <- Y[, 1:k]
  Xk <- X[, 1:k]

  PX <- X %*% l
  PY <- Y %*% l

  Xc <- X %*% sqrt(l)
  Yc <- Y %*% sqrt(l)

  PXk <- PX[, 1:k]
  PYk <- PY[, 1:k]

  Xck <- Xc[, 1:k]
  Yck <- Yc[, 1:k]

  if (exists('SP')) {
    p <- cat[2]
    SPOrig <- SP
    top <- apply(Fc, 1, max)
    SPSize <- dim(SP)
    top <- top * SPOrig
    SPmax <- top * SPOrig
    SPoriginv <- 1 - SPOrig
    Mean <- (sum(Fc[1,]) / p) * SPoriginv
    SPm <- Mean + SPmax
    sumSPm <- apply(SPm, 2, sum)
    diagSum <- diag(sumSPm)
    SP <- as.matrix(SPm) %*% solve(diagSum)


    if (RorC == 'C') {
      Dsp <- diag(apply(SPOrig, 2, sum)) / cat[1]
      SP <- sqrt(solve(Dsp)) %*% as.matrix((t(SP))) %*% Yk
    } else {
      SP <- (t(SP)) %*% PXk
    }
  }


  # Plotting ----
  plot(PXk, main = 'Analysis Fc: Principal coordinates Objects', col = 'blue', pch = 19, cex = 0.5)
  p <- cat[2] / 2
  colors <- rainbow(p)
  for (i in 1:p){
    segments(PXk[i, 1], PXk[i, 2], PXk[p + i, 1], PXk[p + i, 2], col = colors[i])
    text(PXk[i, 1], PXk[i, 2], labels = labelx[i, 1], cex = 0.6, adj = 0)
    text(PXk[i + p, 1], PXk[i + p, 2], labels = labelx[i + p, 1], cex = 0.6, adj = 0)
  }

  if (exists('SP')) {
    if (RorC == 'C') {
      points(SP[, 1], SP[, 2], col = 'red', pch = 16, cex = 0.6)
      for (i in 1:dim(SP)[1]) {
        text(SP[i, 1], SP[i, 2], labels = LabelsSup[i, 1], cex = 0.6, adj = 0)
      }
    }
  } else {
    plot(Xk, main = 'Analysis Fc: Standard coordinates Objects', col = 'blue', pch = 19, cex = 0.5)
    for (i in 1:p){
      segments(Xk[i, 1], Xk[i, 2], Xk[p + i, 1], Xk[p + i, 2], col = colors[i])
      text(Xk[i, 1], Xk[i, 2], labels = labelx[i, 1], cex = 0.6, adj = 0)
      text(Xk[i + p, 1], Xk[i + p, 2], labels = labelx[i + p, 1], cex = 0.6, adj = 0)
    }
    points(SP[, 1], SP[, 2], col = 'red', pch = 16, cex = 0.6)
    for (i in 1:dim(SP)[1]) {
      text(SP[i, 1], SP[i, 2], labels = LabelsSup[i, 1], cex = 0.6, adj = 0)
    }
  }

}

CAFc(C, 2, Zsupsred, Lsupsred, 'C')
