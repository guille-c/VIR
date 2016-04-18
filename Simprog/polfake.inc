      INTEGER SIZE, IDIM, CDIM
      PARAMETER (SIZE=1024)
      PARAMETER (IDIM=SIZE*SIZE)
      PARAMETER (CDIM=(SIZE/2+1)*SIZE)
C Flags are TRUE if the corresponding image is available
      LOGICAL DOCMB, DOFG
C Dimensions of images(s)
      INTEGER MDIM(2)
C Pixel size (degrees)
      DOUBLE PRECISION PIXL(2)
C Reference pixel (map center)
      DOUBLE PRECISION REFP(2)
C Arrays to hold CMB images
      REAL MAPARR(IDIM), QMAP(IDIM), UMAP(IDIM)
C Array to hold foreground image
      REAL FGARR(IDIM)
C Arrays for doing the FFT (image and uv plane)
      REAL MAPWT(IDIM), MAPWTQ(IDIM), MAPWTU(IDIM)
      COMPLEX MAPFT(CDIM), MAPFTQ(CDIM), MAPFTU(IDIM)
      EQUIVALENCE (MAPWT, MAPFT)
      EQUIVALENCE (MAPWTQ, MAPFTQ)
      EQUIVALENCE (MAPWTU, MAPFTU)
C All stored in common
      COMMON /MAP2/ DOCMB, DOFG, MDIM, PIXL, REFP
      COMMON /MAP1/ MAPARR, QMAP, UMAP, FGARR, MAPWT, MAPWTQ, MAPWTU

