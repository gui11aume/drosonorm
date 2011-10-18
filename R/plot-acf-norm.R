plot.acf.norm <- function(core.name, name, intensity, norm, array1.name,
    array2.name="", maxlag=200, cex, graph) {

  out.file <- paste(core.name,"_ACF_", .dtag(), ".png",sep="");
  n.arrays <- ifelse(array2.name == "", 1, 2);

  if (n.arrays == 1) {
    canvas.width <- 256;
    mfrow <- c(1,1);
  }
  else {
    canvas.width <- 1024;
    mfrow <- c(1,4);
  }

  if (graph == "x11") {
    png(
        file = out.file,
        width = canvas.width,
        height = 400);
  }
  if (graph == "ps") {
    bitmap(
        file = out.file,
        width = canvas.width,
        height = 400,
        units = "px",
        point = 1,
        taa = 4);
  }


  ###############################
  ##      "DRY" functions      ##
  ###############################

  acf1array <- function(x, main, col) {
    acf(x=x, lag.max=maxlag, na.action=na.pass,
        ci=0, main=main, ylim=c(0,1), col=col);
  }

  addlines <- function() {
    abline(h=c(.2, .4, .6), lty=2, col="grey");
  }

  plotacf <- function(acf.obj, col) {
    lines(acf.obj$lag[1:maxlag], acf.obj$acf[1:maxlag], col=col);
  }

  print.name <- function() {
    txt <- paste(name, sub("([^%])$", "\\1%", intensity));
    mtext(text=txt, line=0.4, col="grey50", cex=1.5*cex);
  }


  par(mfrow=mfrow, cex=cex);

  if (n.arrays == 1) {
    acf.1 <- acf1array(
        x = norm$M.norm,
        main = paste("ACF of array", array1.name),
        col = "red"
    );

    addlines();
    # Print name and intensity.
    print.name();

    # Print version control.
    .print.vcontrol("right", cex=cex);
  }

  if (n.arrays == 2) {
    # Panel 1.
    acf.1 <- acf1array(
        x = norm$M1.norm,
        main = paste("ACF of array", array1.name),
        col = "red"
    );

    addlines();
    # Print name and intensity.
    print.name();

    # Panel 2.
    acf.2 <- acf1array(
        x = norm$M2.norm,
        main = paste("ACF of array", array2.name),
        col = "green"
    );

    addlines();

    # Panel 3.
    acf.av <- acf1array(
        x = norm$M.norm,
        main = "ACF of mean",
        col="black"
    );

    addlines();

    # Panel 4.
    plot(
        x = c(0, maxlag),
        y = c(0,1),
        type = "n",
        main = "Overlay",
        xlab = "Lag",
        ylab = "ACF"
    );

    plotacf(acf.1, col="red");
    plotacf(acf.2, col="green");
    plotacf(acf.av, col="black");

    # Print version control.
    .print.vcontrol("right", cex=cex);

  }

  dev.off();

}
