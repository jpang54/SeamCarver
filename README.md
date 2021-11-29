# SeamCarver

#### This is a content-aware image resizing algorithm that is able to resize pictures without losing the main components within it. 

#### For example, in a picture with three people surfing and one surfer a large distance away from the others, this algorithm can crop out this distance and bring them closer together while ensuring the picture still looks natural. The algorithm chooses which seams to remove in the picture by calculating the energy of each pixel, where pixels with rapid color gradients have higher energy and are thus less likely to be included in a seam to be removed. This algorithm uses an alternate version of Dijkstra’s algorithm to calculate the lowest energy seam to be removed, optimized to be faster by caching already calculated energy values and to use less memory by representing the digraph implicitly using a 2d array. 
