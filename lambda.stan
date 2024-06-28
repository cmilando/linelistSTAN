  // Function to compute lambda without using rev
  vector lambda(vector curve, vector si) {
    int k = num_elements(si);
    int n = num_elements(curve);
    vector[n - 1] result;

    for (int i = 1; i <= n - 1; ++i) {
      if (i <= k) {
        vector[i] c = head(curve, i);
        vector[i] s = head(si, i);
        for (int j = 1; j <= i; ++j) {
          result[i] += c[j] * s[i - j + 1];
        }
      } else {
        vector[k] c = segment(curve, i - k, k);
        vector[k] s = si;
        for (int j = 1; j <= k; ++j) {
          result[i] += c[j] * s[k - j + 1];
        }
      }
    }
    return result;
  }


