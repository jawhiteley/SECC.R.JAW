##################################################
# Schefferville Experiment on Climate Change (SEC-C)
# R v2.12
##################################################
## tests of outputs to different file formats

txt <- "Hello World ±\n123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 1234567890"
txt <- paste(txt, width = 24)

x1 <- quote( mu*"mol" %.% m^-2 )
x2 <- eval("d^-1")
x3 <- expression(chi^2)
x.text <- "Label Text"

txt2 <- bquote( .(x.text) * " (" * .(x1) %.% .(x2) %*% quote(.(x3)) * ")" )
txt3 <- bquote( paste( .(x.text), " (", 
					  .( bquote(.(x1) %.% .(x2) %*% quote(.(x3)) ) ),
					  ")" 
					  )
			  )

## pdf( file = "~/Desktop/test.pdf" )

plot.new()
# plot.window(c(0, 100), c(0, 100))
mtext( txt3 ,
       adj = 0, side = 3, line = 0, 
       family = "mono", cex = 1, col = "black"
      )
plot( rep(1, 25), 1:25, pch = 1:25 )

## dev.off()
