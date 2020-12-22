# Mean_variance plots

# objects from the _tweedie scripts
meanvar.plot(elements_only)
shape_only
elements_and_shape

options(scipen=99)

par(mfrow=c(1,3))
meanvar.plot(elements_only, mfrow =NULL, mfcol = NULL)
meanvar.plot(shape_only, mfrow =NULL, mfcol = NULL, xlab="Mean (log scale)", yaxt='n')
axis(side=2, at=c(0.001,0.01,0.1), labels=c(0.001,0.01,0.1))
meanvar.plot(elements_and_shape, mfrow =NULL, mfcol = NULL)
