
  // * x = repporting delay (integer)
  // * onset = onset day (integer)
  vector prop(vector x, vector onset, int maxdelay, int cd) {
    int n = num_elements(x);
    vector[n] x1;
    int dem;
    vector[maxdelay] p1;
    vector[maxdelay] p;
    vector[maxdelay] p2;
    vector[maxdelay] result;
    int count = 1;

    // Create logical vector
    for (i in 1:n) {
      if (x[i] <= maxdelay && onset[i] >= cd) {
        x1[count] = x[i];
        count += 1;
      }
    }

    // Truncate x1 to the actual size of valid elements
    vector[count-1] x1_truncated;
    for (i in 1:(count-1)) {
      x1_truncated[i] = x1[i];
    }

    dem = count - 1;

    // Initialize p1
    for (i in 1:maxdelay) {
      p1[i] = 0;
    }

    // Populate p1
    for (i in 1:maxdelay) {
      for (j in 1:dem) {
        if (x1_truncated[j] == (maxdelay - i + 1)) {
          p1[i] += 1;
        }
      }
    }

    // Calculate proportion
    p = p1 / dem;

    // Calculate cumulative sum
    p2[1] = p[1];
    for (i in 2:maxdelay) {
      p2[i] = p2[i-1] + p[i];
    }

    // Calculate result
    for (i in 1:maxdelay) {
      result[i] = 1 - p2[i];
    }

    return result;
  }
