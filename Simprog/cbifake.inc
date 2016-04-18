      INTEGER SIZE, IDIM, CDIM
      PARAMETER (SIZE=2048)
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
C Array to hold CMB image
      REAL MAPARR(IDIM)
C Array to hold foreground image
      REAL FGARR(IDIM)
C Array for doing the FFT (image and uv plane)
      REAL MAPWT(IDIM)
      COMPLEX MAPFT(CDIM)
      EQUIVALENCE (MAPWT, MAPFT)
C All stored in common
      COMMON /MAP2/ DOCMB, DOFG, MDIM, PIXL, REFP
      COMMON /MAP1/ MAPARR, FGARR, MAPWT

