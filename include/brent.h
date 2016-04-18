#ifdef __cplusplus
extern "C" {
#endif  

float brent(float ax, float bx, float cx, float (*f)(float), float tol,
	    float *xmin);

#ifdef __cplusplus
}
#endif  

