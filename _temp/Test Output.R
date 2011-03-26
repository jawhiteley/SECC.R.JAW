##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# R v2.12.1
##################################################
## tests of outputs to different file formats

txt <- "Hello World ±\n123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 1234567890"

pdf( file = "~/Desktop/test.pdf" )

plot.new()
# plot.window(c(0, 100), c(0, 100))
mtext( paste(txt, width = 24), 
       adj = 0, side = 3, line = 0, 
       family = "mono", cex = 1, col = "black"
      )
plot( rep(1, 25), 1:25, pch = 1:25 )

dev.off()
