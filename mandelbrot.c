#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

// picture resolution
static const int ImageWidth  = 1000;
static const int ImageHeight = 1000;
static const double CxMin    = -2.5;
static const double CxMax    = 1.5;
static const double CyMin    = -2.0;
static const double CyMax    = 2.0;

int main()
{
    double PixelWidth=(CxMax-CxMin)/ImageWidth; /* scaled x coordinate of pixel (must be scaled to lie somewhere in the Mandelbrot X scale (-2.5, 1.5) */
    double PixelHeight=(CyMax-CyMin)/ImageHeight;/* scaled y coordinate of pixel (must be scaled to lie somewhere in the Mandelbrot Y scale (-2.0, 2.0) */

    const int MaxColorComponentValue=1<<8;
    typedef unsigned char pixel_t[3]; // colors [R, G ,B]
    pixel_t *pixels = malloc(sizeof(pixel_t)*ImageHeight*ImageWidth);

    //write the header to the picture file
    FILE* fp;
    #pragma omp parallel shared(pixels)
    {
        int iY;
        #pragma omp for
        for(iY=0; iY<ImageHeight; iY++) {
            double Cy=CyMin + iY*PixelHeight;
            if (fabs(Cy)< PixelHeight/2) {
                Cy=0.0; // Main antenna
            }
            int iX;
            for(iX=0; iX<ImageWidth; iX++) {
                double Cx=CxMin + iX*PixelWidth;
                double Zx=0.0;
                double Zy=0.0;
                double Zx2=Zx*Zx;
                double Zy2=Zy*Zy;
                /* */
                int Iteration;
                const int IterationMax=150;
                const double Bailout=2; // bail-out value
                const double Circle_Radius=Bailout*Bailout; // circle radius

                for (Iteration=0; Iteration<IterationMax && ((Zx2+Zy2)<Circle_Radius); Iteration++) { // 
                    Zy=2*Zx*Zy + Cy;
                    Zx=Zx2-Zy2 + Cx;
                    Zx2=Zx*Zx;
                    Zy2=Zy*Zy;
                };

                if (Iteration==IterationMax) {
                    //  interior of Mandelbrot set = black
                    pixels[iY*ImageWidth + iX][0] = 0;
                    pixels[iY*ImageWidth + iX][1] = 0;
                    pixels[iY*ImageWidth + iX][2] = 0;
                }
                //
                else {
                    pixels[iY*ImageWidth + iX][0] = ((double)(Iteration-log2(log2(sqrt(Zx2+Zy2))))/IterationMax) * MaxColorComponentValue;
                    pixels[iY*ImageWidth + iX][1] = 0;
                    pixels[iY*ImageWidth + iX][2] = 0;
                }
            }
        }
    }
    //return 0; // TO BENCHMARK without free()/file IO
    //write the header to the picture file
    fp= fopen("MandelbrotSet.ppm","wb");
    fprintf(fp,"P6\n %s\n %d\n %d\n %d\n","# no comment",ImageWidth,ImageHeight,MaxColorComponentValue);
    for(int iY=0; iY<ImageHeight; iY++)
        for(int iX=0; iX<ImageWidth; iX++)
            fwrite(pixels[iY*ImageWidth + iX],1,sizeof(pixel_t),fp);
    fclose(fp);
    free(pixels);
    //  stop_timer ( );
    //
    //  printf("Elapsed time: %lf\n",elapsed_time ( ));
    return 0;
}
