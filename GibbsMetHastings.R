
GenerateDataFromDistr <- function(a, b, c, tau = 1) {
  # Generation of data from the NHPP distribution.
  # 
  # Args:
  #   a, b, c: Parameters in distribution.
  #   tau: End time.
  #   
  # Returns:
  #   A vector of samples from model
  tauLambda = LargeLambdaEstimate(tau, a, b, c)
  numpoints = rpois(1, tauLambda)
  # vector to be returned
  x = replicate(numpoints, GenerateOneSamplePoisson(tauLambda, a, b, c))
  return(x)
}

GenerateOneSamplePoisson <- function(tauLambda, a, b, c) {
  # Finds x for the cumulative distribution equal a probability.
  # 
  # Args:
  #   tauLambda: Large lambda estimate for tau.
  #   
  # Returns:
  #   A sample from poisson distribution.
  
  # Cumulative probability.
  u = runif(1) 
  # Starting time.
  x = 0
  tolerance = 0.0001
  stepSize = 0.1
  largeLambdaValue = 0
  direction = 1
  
  while (abs(largeLambdaValue - u) > tolerance)  {
    x = x + direction*stepSize
    # Check that x is still above 0
    if (x <= 0) {
      direction = 1
      x = stepSize
    }
    largeLambdaValue = LargeLambdaEstimate(x, a, b, c)/tauLambda
    
    # Going left and pass the point
    if ((u > largeLambdaValue) && (direction == -1)) {
      stepSize = stepSize/2
      direction = 1
    }
    
    # Going right and pass the point
    if ((u < largeLambdaValue) && (direction == 1)) {
      stepSize = stepSize/2
      direction = -1
    }
  }
  return(x)
}

GenerateDataNumInverse <- function(a, b, c, NUMPOINTS) {
  # Generation of data from a NHPP model.
  # Done by finding inverse numerical.
  #
  # Args:
  #   a, b, c: Parameters in model.
  #
  # Returns:
  #   A vector of samples from model.
  
  # One too many to normalize the sample between 0 and 1.
  x = rep(0, NUMPOINTS + 1)
  for (i in 1:length(x)) {
    x[i] = GeneratOneSampleNumInverse(a, b, c)
  }
  x = sort(x)
  # Normalizing times between 0 and 1.
  x = x/x[length(x)]
  return(x[1:NUMPOINTS])
}

GeneratOneSampleNumInverse <- function(a, b, c) {
  # Generate t value for a value u. Numerically.
  #
  # Args:
  #   a, b, c: Parameters in NHPP.
  #
  # Returns:
  #   One sample from NHPP
  
  # Cumulative probability
  u = runif(1)
  
  # Starting time.
  x = 0
  tolerance = 0.0001
  stepSize = 0.1
  largeLambdaValue = 0
  direction = 1
  
  while (abs(largeLambdaValue - u) > tolerance)  {
    x = x + direction*stepSize
    # Check that x is still above 0
    if (x <= 0) {
      direction = 1
      x = stepSize
    }
    largeLambdaValue = LargeLambdaEstimate(x, a, b, c)
    
    # Going left and pass the point
    if ((u > largeLambdaValue) && (direction == -1)) {
      stepSize = stepSize/2
      direction = 1
    }
    
    # Going right and pass the point
    if ((u < largeLambdaValue) && (direction == 1)) {
      stepSize = stepSize/2
      direction = -1
    }
  }
  
  return(x)
}

ScrambleData <- function(xSampleInit, statS1, statS2,
                         NUMPOINTS, nIterations) {
  # Scramble data by picking random x, but keeping the sufficient statistics.
  #
  # Args: 
  #   xSampleInit: Data to scramble.
  #   statS1: Sum of data.
  #   statS2: Prod of data.
  #   NUMPOINT: Number of points in data.
  #   nIterations: Number of iterations for scramble.
  #   
  # Returns:
  #   A vector of scrambled data.

  myData = xSampleInit
  for(i in 1:nIterations){
    xThreeUnique = DrawThreeUniqueX(NUMPOINTS)
    # Calculate sum and prod for these x's.
    sumX = myData[xThreeUnique[1]] + myData[xThreeUnique[2]] + myData[xThreeUnique[3]]
    prodX = myData[xThreeUnique[1]] * myData[xThreeUnique[2]] * myData[xThreeUnique[3]]
    
    newX3 = runif(1)
    if (sumX < 1) {
      newX3 = newX3*sumX
    }
    if (IsValid(newX3, sumX, prodX)){
      roots = FindRoots(newX3, sumX, prodX)
      if (CheckZeroToOne(roots[1]) && CheckZeroToOne(roots[2])) {
        # Update myData
        myData[xThreeUnique[1]] = roots[1]
        myData[xThreeUnique[2]] = roots[2]
        myData[xThreeUnique[3]] = newX3 
      }
    }
  }
  return(myData)
}

SamplerGibbsMet <- function(xSampleInit, statS1, statS2,
                            nIterations, convergencePlot = FALSE) {
  # A Gibbs sampler with one metropolis hastings step.
  #
  # Args:
  #   xSampleInit: Starting vector for the sample.
  #   statS1: Sum of data.
  #   statS2: Prod of data.
  #   nIterations: number of iterations before returning sample.
  #   convergerncePlot: Check for convergence.
  # 
  # Returns:
  #   Sample given statS1 and statS2. A vector.
  
  mySample = xSampleInit
  
  yConvPlotValues = 0
  if (convergencePlot) {
    yConvPlotValues = rep(0, nIterations)
  }
  
  for (i in 1:nIterations) {
    xThreeUnique = DrawThreeUniqueX(length(xSampleInit))
    # Calculate sum and prod for these x's.
    sumX = mySample[xThreeUnique[1]] + mySample[xThreeUnique[2]] + mySample[xThreeUnique[3]]
    prodX = mySample[xThreeUnique[1]] * mySample[xThreeUnique[2]] * mySample[xThreeUnique[3]]
    
    x3Current = mySample[xThreeUnique[3]]
    # Draw proposal X3
    x3 = runif(1)
    # Check validity of x3
    if (sumX < 1) {
      x3 = x3*sumX
    }
    if (IsValid(x3, sumX, prodX)) {
      # Calculate x1 and x2
      roots = FindRoots(x3, sumX, prodX)
      # Check roots for legal values.
      if(CheckZeroToOne(roots[1]) && CheckZeroToOne(roots[2]) && !is.nan(roots[1]) && !is.nan(roots[2])) {
        # Calculate acceptance
        alpha = min(1, Distribution(x3, sumX, prodX)/Distribution(x3Current, sumX, prodX))
        # Acceptance step
        u = runif(1)
        if(u <= alpha) {
          # Update mySample
          mySample[xThreeUnique[1]] = roots[1]
          mySample[xThreeUnique[2]] = roots[2]
          mySample[xThreeUnique[3]] = x3 
#           if (convergencePlot) {
#             yConvPlotValues[i] = yConvPlotValues[i] + 1
#           }
        }
      }
    }
    
    if (convergencePlot) {
      yConvPlotValues[i] = median(mySample)
    }
  }
  if (convergencePlot) {
    print(sum(yConvPlotValues)/nIterations)
    plot(seq(1, nIterations, 1), yConvPlotValues, type="l")
  }
  return(mySample)
}

FindXByInverseF <- function(sumX, prodX) {
  # Generating a value x by finding the inverse of F(x).
  # 
  # Args:
  # 
  # Returns:
  #   A value for the inverse F.
  
  # Cumulative probability
  u = runif(1)
  
  # Starting time.
  x = 0
  tolerance = 0.0001
  stepSize = 0.1
  integralValue = 0
  direction = 1
  
  while (abs(integralValue - u) > tolerance)  {
    x = x + direction*stepSize
    # Check that x is still above 0
    if (x <= 0) {
      direction = 1
      x = stepSize
    }
    integralValue = CalcCumulativeDistr(x, sumX, prodX)
    
    # Going left and pass the point
    if ((u > integralValue) && (direction == -1)) {
      stepSize = stepSize/2
      direction = 1
    }
    
    # Going right and pass the point
    if ((u < integralValue) && (direction == 1)) {
      stepSize = stepSize/2
      direction = -1
    }
  }
  return(x)
}

CalcCumulativeDistr <- function(x, sumX, prodX) {
  # Using the trapez method to calculate cumulative functon for given x.
  #
  # Args:
  #   x: End point for the integral.
  #   sumX: Sum of the x's
  #   prodX: Prod of the x's
  # 
  # Returns:
  #   value for integral.

  # Number of steps.
  n = 100
  # Start cannot be 0. Trapez function becomes divide by zero for some b = 0. 
  start = 0.0001
  end = x
  stepSize = (end-start)/n
  x = seq(start, end, stepSize)
  y = rep(0, length(x))
  y = 1/(x*sqrt((sumX - x)^2 - 4*prodX/x))
  tmp = (y > 0)
  y = y*tmp
  return((stepSize/2)*(y[1] + y[n+1] + 2*sum(y[2:n])))  
  
}

DrawThreeUniqueX <- function(numPoints) {
  # Draw three unique values from 1 to NUMPOINTS.
  #
  # Args:
  #   NUMPOINTS: End value for draw intervall.
  #
  # Retuns:
  #   A vector of length 3 with unique values.
  
  return(sample(numPoints, 3))
}

IsValid <- function(value, a, b) {
  # Check to see if x1 is a valid value for our density function.
  #
  # Args:
  #   value: Value to check for validity.
  #   a: Sum.
  #   b: Product.
  #
  # Returns:
  #   A boolean value 1 or 0.
  
  return(((value*((a - value)^2) - 4*b) > 0))
}

FindRoots <- function(x3, a, b) {
  # Finding roots of second order function.
  #
  # Args:
  #   x3: Proposal value.
  #   a: Sum of the three x's.
  #   b: Product of the three x's.
  #
  # Returns:
  #   A vector of length 2 containing the two roots.
  
  roots = rep(0, 2)
  roots[1] = ((a - x3) + sqrt((a - x3)^2 - 4*b/x3))/2 
  roots[2] = ((a - x3) - sqrt((a - x3)^2 - 4*b/x3))/2
  return(roots)
}

CheckZeroToOne <- function(value) {
  # Check to see if a value is between 0 and 1.
  #
  # Args:
  #   value: Value to check.
  #
  # Returns:
  #   A boolean value 1 or 0.
  return((value <= 1) && (value >= 0))
}

Distribution <- function(x3, a, b) {
  # Calculating probability for a value x3.
  #
  # Args:
  #   x3: Variable.
  #   a: Sum of the three x's
  #   b: Prod of the three x's
  #
  # Returns:
  #   A probability value.
  return(1/(x3*sqrt(((a-x3)^2) - 4*b/x3)))
} 

LogLikelihoodNhpp <- function(x) {
  # Calculating value of log likelihood to NHPP.
  #
  # Args:
  #   x: Vector of parameters in log likelihood.
  #
  # Returns:
  #   A value for the log likelihood.
  
  a = x[1]
  b = x[2]
  c = x[3]
  dataPoints = length(data)
  # This the negative log likelihood. This is for finding max(min).
  L = -dataPoints*log(a) - dataPoints*log(b) - (b-1)*sum(log(data)) - c*sum(data) + a*b*TrapezIntegral(1, b, c)
  return(L)
}

TrapezIntegral <- function(end, b, c) {
  # Trapez method to solve an integral.
  #
  # Args: 
  #   b: Value of parameter.
  #   c: Value of parameter.
  #
  # Returns:
  #   Value of the integral.
  
  # Step size
  n = 100
  # Start cannot be 0. Trapez function becomes divide by zero for some b = 0. 
  start = 0.0001
  stepSize = (end-start)/n
  x = seq(start, end, stepSize)
  y = TrapezFunction(x, b, c)
#   print(y)
#   if (y[1] == Inf) {
#     print("start")
#     print(x[1])
#     print(b)
#     print(c)
#     print("slutt")
#   }
  return((stepSize/2)*(y[1] + y[n+1] + 2*sum(y[2:n])))
}

TrapezFunction <- function(x, b, c) {
  # Function to itegrate over.
  #
  # Args:
  #   x: Data point.
  #   b: Value for b parameter.
  #   c: value for c parameter.
  #
  # Returns:
  #   A value for given data point.
  return((x^(b-1))*exp(c*x))
}

SmallLambda <- function(x, a, b, c) {
  # Returns value of lambda function for given parameters
  # 
  # Args:
  #   x: Time value.
  #   a: Parameter value.
  #   b: Parameter value.
  #   c: Parameter value.
  return(a*b*(x^(b-1))*exp(c*x))
}

LargeLambdaEstimate <- function(x, a, b, c) {
  # Estimate of the large lambda function
  #
  # Args:
  #   x: Time value. Not a vector.
  #   a: Estimated value for parameter.
  #   b: Estimated value for parameter.
  #   c: Estimated value for parameter.
  #
  # Returns:
  #   A value for the lambda estimate.
  return(a*b*TrapezIntegral(x, b, c))
}

GreenwoodStatistic <- function(x) {
  # Calculating Greenwood statistic value for given data set. Two sided.
  #
  # Args:
  #   x: Given data set of time.
  #
  # Returns:
  #   Greenwood value
  
  x = sort(x)
  # Transformed times
  v = rep(0, length(x)+1)
  
  greenwoodSum = 0
  
  # Large lambda estimate for xn
  largeLambdaTau = LargeLambdaEstimate(1, aEstimated, bEstimated, cEstimated)
  # Calculate tranformed times
  for (i in 1:length(x)) {
    largeLambdaI = LargeLambdaEstimate(x[i], aEstimated, bEstimated, cEstimated)
    v[i+1] = largeLambdaI/largeLambdaTau
    greenwoodSum = greenwoodSum + (v[i+1] - v[i])^2
  }
  greenwoodSum = greenwoodSum + (1 - v[length(v)])^2
  return(greenwoodSum)
}

CalcStatistics <- function(x) {
  # Calculation of tranformed times and different statistics.
  # 
  # Args:
  #   x: Data given as time.
  
  x = sort(x)
  # Tranformed times
  v = rep(0, length(x))
  largeLambdaTau = LargeLambdaEstimate(1, aEstimated, bEstimated, cEstimated)
  laplaceSum = 0
  cramerSum = 0
  dPlus = -100
  dMinus = -100
  for(i in 1:length(x)) {
    largeLambdaI = LargeLambdaEstimate(x[i], aEstimated, bEstimated, cEstimated)
    v[i] = largeLambdaI/largeLambdaTau
    laplaceSum = laplaceSum + v[i] - 1/2
    cramerSum = cramerSum + (v[i] - ((2*i - 1)/(2*length(x))))^2
    newDplus = (i/length(x) - v[i])
    if (newDplus > dPlus) {
      dPlus = newDplus
    }
    newDMinus = v[i] - ((i-1)/length(x))
    if (newDMinus > dMinus) {
      dMinus = newDMinus
    }
  }
  laplaceSum = sqrt(12/length(x))*laplaceSum
  cramerSum = cramerSum + 1/(12*length(x))
  return(c(cramerSum, laplaceSum, max(dPlus, dMinus)))
}


FindLogLikeliHoodEstimators <- function() {
  # Finding estimators from log likelihood
  #
  # Args:
  #
  # Returns:
  #   A list of parameters.
  
  # Solving numerically.
  result = optim(c(1,1,1), LogLikelihoodNhpp)
  return(result$par)
}

# Test functions
TestCheckZeroToOne <- function() {
  # Check to see if CheckZeroToOne behaves correctly.
  
  print(CheckZeroToOne(-6) == 0)
  print(CheckZeroToOne(6) == 0)
  print(CheckZeroToOne(0.5) == 1)
}

TestSums <- function() {
  # To check that two sums are equal. Two ways of finding x1 + x2.
  
  # Environment setup
  n = 20
  x = runif(n)
  aMark = sum(x[1:3])
  x3 = runif(1)
  a = sum(x) - x3 - sum(x[4:n])
  print(a == (aMark - x3))  
}

# Test runs
TestCheckZeroToOne()
TestSums()

# Running of script

# Iterations between samples. 
NITERATIONS = 1500
# Number of samples
NSAMP = 100000
# Data variable for log likelihood
data = 0
# Model parameters
aModel = 160
bModel = 2
cModel = -3
tau = 1
# Estimated parameters. 
# Estimated parameters will be the same for all samples because sum and product is constant.
aEstimated = 0
bEstimated = 0
cEstimated = 0

# Generate data
xData = GenerateDataFromDistr(aModel, bModel, cModel, tau) 
NUMPOINTS = length(xData)
data = xData
t = seq(from = 0, to = 1, by = 0.01)
lambdaValues = SmallLambda(t, aModel, bModel, cModel)
plot(t, lambdaValues, type="l", main="", ylab = "Lambda values", xlab = "Time values", cex.lab=1.5)
hist(xData, breaks=5, main="", xlab="Time values", cex.lab=1.5)
param = FindLogLikeliHoodEstimators()
aEstimated = param[1]
bEstimated = param[2]
cEstimated = param[3]
statS1 = sum(xData)
statS2 = prod(xData)
greenwoodObs = GreenwoodStatistic(xData)
statisticsObs = CalcStatistics(xData)
cramerObs = statisticsObs[1]
laplaceObs = statisticsObs[2]
kolmogorovObs = statisticsObs[3]


# Used for p-value
greenwoodBelow = 0
greenwoodAbove = 0
laplaceBelow = 0
laplaceAbove = 0
cramerBelow = 0
cramerAbove = 0
kolmogorovAbove = 0


# In case of burn in
burnIterations = 4000
xSample = SamplerGibbsMet(xData, statS1, statS2, burnIterations)

# Statistics storage
greenwoodSamples = rep(0, NSAMP)
laplaceSamples = rep(0, NSAMP)
cramerSamples = rep(0, NSAMP)
kolmogorovSamples = rep(0, NSAMP)


# Simulations
for (i in 1:NSAMP) {
  # Generate sample
  xSample = SamplerGibbsMet(xSample, statS1, statS2, NITERATIONS)
  sampleStatS1 = sum(xSample)
  sampleStatS2 = prod(xSample)
  
  ## Check sufficient statistics for correct values.
  #validStatS1 = (abs(sampleStatS1 - statS1) < 0.1) 
  #validStatS2 = (abs(sampleStatS2 - statS2) < 0.1)
  
  # Calculate greenwood statistic
  greenwood = GreenwoodStatistic(xSample)
  greenwoodSamples[i] = greenwood 
  
  if (greenwood >= greenwoodObs) {
    greenwoodAbove = greenwoodAbove + 1
  }
  
  if (greenwood <= greenwoodObs) {
    greenwoodBelow = greenwoodBelow + 1
  }
  
  statistics = CalcStatistics(xSample)
  
  # Laplace
  laplace = statistics[2]
  laplaceSamples[i] = laplace
  
  if (laplace >= laplaceObs) {
    laplaceAbove = laplaceAbove + 1
  }
  
  if (laplace <= laplaceObs) {
    laplaceBelow = laplaceBelow + 1
  }
  
  # Cramer
  cramer = statistics[1]
  cramerSamples[i] = cramer
  
  if (cramer >= cramerObs) {
    cramerAbove = cramerAbove + 1
  }
  
  # Kolmogorov
  kolmogorov = statistics[3]
  kolmogorovSamples[i] = kolmogorov
  
  if (kolmogorov >= kolmogorovObs) {
    kolmogorovAbove = kolmogorovAbove + 1
  }
  
}

# Plot of a sample
hist(xSample, breaks = 5, main="", xlab="Time values", cex.lab=1.5)

# Calculate p-value
pValueGreenwood = 2*min(greenwoodAbove/NSAMP, greenwoodBelow/NSAMP)
pValueLaplace = 2*min(laplaceAbove/NSAMP, laplaceBelow/NSAMP)
pValueCramer = cramerAbove/NSAMP
pValueKolmogorov = kolmogorovAbove/NSAMP
# Plot of greenwood statistics
breaks = 200
hist(greenwoodSamples, breaks=breaks, main="", xlab="Greenwood values", cex.lab=1.5)
abline(v = greenwoodObs, col="red")
hist(laplaceSamples,breaks=breaks, main="", xlab="Laplace values", cex.lab=1.5)
abline(v = laplaceObs, col="red")
hist(cramerSamples,breaks=breaks, main="", xlab="Cramer values", cex.lab=1.5)
abline(v = cramerObs, col="red")
hist(kolmogorovSamples, breaks=breaks, main="", xlab="Kolmogorov values",  cex.lab=1.5)
abline(v = kolmogorovObs, col="red")

# Save image
save.image(file="myData.RData")

print("Done")
