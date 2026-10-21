# Create the nimble function object
cumulative_trapz <- nimbleFunction(
  # run = function(f_values = double(1), dx = double(1)) { ... }
  # NOTE: Assuming dx is a vector of segment widths, where length(dx) == length(f_values) - 1.
  # If dx is a scalar representing constant spacing, see the alternative below.
  run = function(f_values = double(1), dx = double(0)) {
    n_segments <- length(f_values) - 1
    # Check for correct input lengths
    
    # 1. Initialize the result vector (length N)
    # The first element of the cumulative integral is always 0 (integral up to the start)
    cumulative_integral <- numeric(n_segments + 1)
    cumulative_integral[1] <- 0
    
    # 2. Loop through each segment to calculate the area and accumulate it
    for (i in 1:n_segments) {
      # Trapezoid area for segment 'i': 0.5 * (height1 + height2) * width
      # Height 1 is f_values[i], Height 2 is f_values[i+1]
      # Width is dx_values[i]
      trapezoid_area <- 0.5 * (f_values[i] + f_values[i+1]) * dx
      
      # Accumulate the area: new value = previous value + current segment area
      cumulative_integral[i+1] <- cumulative_integral[i] + trapezoid_area
    }
    
    # Return the vector of cumulative integral values
    return(cumulative_integral)
    returnType(double(1))
  }
)
