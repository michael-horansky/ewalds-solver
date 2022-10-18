# ewalds-solver
Brute-force solver that finds scattered wave-vectors for a crystal lattice and an incident wave

Requirements: Python 3, numpy

Notes:
  -If the lattice vectors have a factor of a, you can ignore it, as long as the incident wave-vector has a factor of 1/a (and the scattered wave-vectors then gain the same factor).
  -The equation 2 k_in G + G^2 = 0 is being checked as | 2 k_in G + G^2 | < zero_threshold, where zero_threshold has been set to 1e-06. Feel free to adjust this value to your needs.
  -The trivial solution k_out = k_in is omitted for the sake of brevity.
  -O(n) = n^3 where n is proportional to |k_in|.|a_1|. To minimize computation time, input small crystals and small waves.

hello world
