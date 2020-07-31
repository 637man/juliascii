
// * * *   J U L I A S C I I   F R A C T A L   * * *
// Written by 637man originally in 2018/5, this cleaner version 2020/5
// Compile with just: gcc juliascii.c -Wall -O3 -lpthread -o juliascii

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h> // Always use threads, makes no sense to be singlethreaded.

#define DEFX 0.0f /**< Default center Z real component*/
#define DEFY 0.0f /**< Default center Z imaginary component*/
#define DEFCX -0.8f /**< Default C real component*/
#define DEFCY 0.156f /**< Default C imaginary component*/
#define DEFZOOM 1.0f /**< Default zoom level*/
#define DEFITERS 64 /**< Default amount of iterations*/
#define DEFROTATE 0 /**< Default terminal output rotation mode*/
#define DEFTERMW 80 /**< Default terminal width (40 minimum)*/
#define DEFTERMH 50 /**< Default terminal height (16 minimum)*/
#define DEFPALCH 0 /**< Default ANSI palette*/
#define DEFGRADCH 1 /**< Default color gradient*/
#define DEFEXP 2 /**< Default exponent*/

typedef float FLOAT; // for testing various precisions
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint; // comment for ARM32 (has no effect AFAIK)
//typedef unsigned long uint; // uncomment for ARM32

/**
 * A structure to hold the screen buffer and associated metadata.
 * Supports up to 4 milliard pixels in each dimension on x86_64, 65536 on ARM.
 * Uses 16-bit color depth.
 */
typedef struct {
	uint width; /**< Framebuffer width in pixels*/
	uint height; /**< Framefuffer height in pixels*/
	//uchar rbits; /**< Amount of bits for red component*/
    //uchar gbits; /**< Amount of bits for green component*/
    //uchar bbits; /**< Amount of bits for blue component*/
    ushort** screen; /**< RGB565 format, conversion done by putpixel()*/
} Screen;

/**
 * Constructs a new screen page with all pixels black.
 *
 * @param width     The width of the screen in pixels
 * @param height    The height of the screen in pixels
 * @return          A pointer to the newly constructed Screen structure
 */
Screen* newscreen(uint width, uint height){
	Screen* scr = (Screen*)calloc(1,sizeof(Screen));
	scr->width = width;
	scr->height = height;
	scr->screen = (ushort**)calloc(height,sizeof(ushort*));
	for(int i = 0; i < height; ++i){
		scr->screen[i] = (ushort*)calloc(width,sizeof(ushort));
	}
	return scr;
}

/**
 * Destructs a given Screen structure.
 *
 * @param scr   the screen to be destroyed
 */
void freescreen(Screen* scr){
    for(int i = 0; i < scr->height; ++i){
		free(scr->screen[i]);
	}
    free(scr->screen);
	free(scr);
}

/**
 * Sets a pixel in the screen to the color defined by the last 3 arguments
 * in a RGB8 fashion (unsigned char), converting them to the RGB565 format.
 * Attempting to set a pixel outside of the screen will do nothing.
 *
 * @param scr   The Screen structure pointer where to write to
 * @param x     Column index
 * @param y     Row index
 * @param r     Red component
 * @param g     Green component
 * @param b     Blue component
 */
void setpixel(Screen* scr, uint x, uint y, uchar r, uchar g, uchar b){
	if (y > scr->height-1 || x > scr->width-1) return; // in bounds check
	// TODO make use of src->r/g/bbits?
	ushort rbits = (r << 8) & 0xF800; // 1111 1000 0000 0000
    ushort gbits = (g << 3) & 0x07E0; // 0000 0111 1110 0000
	ushort bbits = b >> 3;
	scr->screen[y][x] = rbits | gbits | bbits;
}

/**
 * Prints a given Screen buffer in an ANSI-art way using the provided palette.
 * Used primarily for debugging, only green component is taken.
 *
 * @param scr   The Screen structure to ANSI-fy
 * @param tw 	The amount of characters per line of the terminal
 * @param th 	The amount of lines of the terminal
 * @param pal   The characters sorted by "luminance"
 * @param plen  The amount of characters in palette
 */
void putscreen(Screen* scr, ushort tw, ushort th, const char* pal, uchar plen){
	ushort printed = 0;
	for(FLOAT y = 0; y < scr->height; y += scr->height/(FLOAT)(th-4)){
		for(FLOAT x = 0; x < scr->width; x += scr->width/(FLOAT)tw){
			putchar(pal[((scr->screen[(uint)y][(uint)x]&0x07E0)>>5)%(plen-1)]);
			++printed;
			if (printed >= tw){
				printed = 0;
				break; // fix for some inaccurracies (namely 132 CPL)
			}
		}
		puts("");
		// not desired under Windows with full tw set, but prevents previous
		// output from messing up after resizing the terminal window
	}
}

/**
 * Same as putscreen() except a filename is taken in addition
 * and the Screen is writen to the file with that name.
 *
 * @param scr   	The Screen structure to ANSI-fy
 * @param fname 	The file name to be written to
 * @param tw 		The amount of characters per line of the terminal
 * @param th 		The amount of lines of the terminal
 * @param pal   	The characters sorted by "luminance"
 * @param plen  	The amount of characters in palette
 */
void savescreen(Screen* scr, const char* fname, ushort tw, ushort th,
		const char* pal, uchar plen){
	FILE* out = fopen(fname, "w");
    for(uint y = 0; y < scr->height; y += scr->height/(FLOAT)th){
		for(uint x = 0; x < scr->width; x += scr->width/(FLOAT)tw){
  		//for(uint y = 0; y < scr->height; y += scr->height/tw){
			fprintf(out, "%c", pal[
					((scr->screen[(uint)y][(uint)x] & 0x07E0) >> 5) % (plen - 1)
				]);
		}
		fprintf(out, "\n");
	}
    fclose(out);
}

/**
 * Wraps drawjulia() arguments, as both functions pthread_create() (POSIX)
 * and thrd_create() can only pass 1 argument to a function. Stupid, ain't it?
 * The construction and destruction of Drawargs is trivial.
 */
typedef struct {
	Screen* scr; /**< Shared framebuffer between threads*/
	int imax; /**< Iteration limit*/
	FLOAT* delims; /**< Individual range for each thread*/
	FLOAT cx; /**< Real added constant component*/
    FLOAT cy; /**< Imaginary added constant component*/
    char exp; /**< Exponent for Z*/
    char grad; /**< Color gradient choice*/
	uint thr; /**< Thread number, denotes which part to draw*/
	uint thrnum; /**< Total threads, denotes how many parts there are*/
} Drawargs;

/**
 * Drawargs setter. There's also thread numerator and denominator.
 * These are used to determine the range within a thread.
 *
 * @param da    The Drawargs structure to configure
 * @param scr	The Screen structure, where to draw the fractal
 * @param imax	Iteration limit
 * @param xmin	Real component of the sequence initial
 * @param ymin	Imaginary component of the sequence initial
 * @param xmax	Real component of the sequence final
 * @param ymax	Imaginary component of the sequence final
 * @param cx	Real component of the constant added
 * @param cy	Imaginary component of the constant added
 * @param exp   Exponent for Z
 * @param grad  Selected gradient
 * @param t     Thread numerator
 * @param tn    Thread denominator
 */
void setdrawargs(Drawargs* da,Screen *scr,uint imax,FLOAT xmin,FLOAT ymin,
FLOAT xmax,FLOAT ymax,FLOAT cx,FLOAT cy,char exp,char grad,uint t,uint tn){
	da->scr = scr;						da->imax = imax;
	da->delims[0] = xmin;				da->delims[1] = ymin;
    da->delims[2] = xmax;				da->delims[3] = ymax;
	da->cx = cx;						da->cy = cy;
	da->exp = exp;						da->grad = grad;
	da->thr = t;                        da->thrnum = tn;
}

/**
 * The core procedure (re)calculating the Julia fractal. Adapted from
 * https://www.root.cz/clanky/fraktaly-v-pocitacove-grafice-x/#k05.
 * Enabled arbitrary integer exponent in the range of char, very yanky results
 * because of odd parity of polynomial functions with odd degree. Still contains
 * some Czech comments.
 *
 * @param argsv  the Drawargs structure address
 */
void drawjulia(void* argsv){ // extract the ordinary arguments
	Drawargs* args = (Drawargs*) argsv;
	Screen* scr = args->scr;			uint imax = args->imax;
	FLOAT xmin = args->delims[0];		FLOAT ymin = args->delims[1];
	FLOAT xmax = args->delims[2];		FLOAT ymax = args->delims[3];
	FLOAT cx = args->cx;				FLOAT cy = args->cy;
	char exp = args->exp;				uchar grad = args->grad;
	uint thr = args->thr;               uint thrnum = args->thrnum;

    FLOAT zx, zy, zx2, zy2, zx0, zy0; // slozsky komplexnii promjennee Z a Z^n
    uint x, y, i, j; // poczitadla sloupcuu a rzsaadkuu v pixmapje a iterace
	FLOAT lim;

    uint ystart = thr;
	uint yend = scr->height;
	uint step = thrnum;

    zy0=ymin;
	lim = 1.0f;
	for(i = 0; i < exp; ++i){
		lim *= 2.0f;
	}

    for (y = ystart; y < yend; y+=step) {   // pro vszechny rzsaadky v pixmapje
        zx0=xmin;
        for (x = 0; x < scr->width; ++x) {  // pro vszechny pixely na rzsaadku
            zx=zx0;                         // nastavit poczaatecznii Z(0)
            zy=zy0;
            for (i = 0; i < imax; ++i) {	// iteracznii smyczka
                zx2=1;
                zy2=1;
				if(exp > 0){                // vyypoczet mocnin slozsek Z
					for(j = 0; j < exp; ++j){
                        zx2*=zx;
                		zy2*=zy;
					}
				} else if (exp < 0){        // negative exponent still exponent
                    for(j = 0; j > -exp; ++j){
                        zx2*=zx;
                		zy2*=zy;
					}
                    zx2=1.0f/zx2;
                	zy2=1.0f/zy2;
				}
                if (zx2+zy2>lim) break;		// przsekroczenii meze divergence
                zy=2.0f*zx*zy+cy;           // vyypoczet Z(n+1)
                zx=zx2-zy2+cx;
            }
			switch(grad){
                case 0: // grayscale
setpixel(scr, x, y, i, i, i);
					break;
				case 1: // blue outside
setpixel(scr, x, y, i * 3, i * 2, i);
					break;
                case 2: // sepia outside
setpixel(scr, x, y, i, i * 2, i * 3);
					break;
                case 3: // green outside
setpixel(scr, x, y, i, i * 2, i);
					break;
				case 4: // exponential blue (psychedelic)
setpixel(scr, x, y, i*i*i, i*i, i);
					break;
                case 5: // exponential red (psychedelic)
setpixel(scr, x, y, i, i*i, i*i*i);
					break;
            	case 6: // subtle noise on red and blue, mostly green
setpixel(scr, x, y, i * 256 / random(), i, i * 256 / random());
					break;
                case 7: // pastel blue-green-red
setpixel(scr, x, y, i+128, i+64, i);
					break;
				case 8: // pastel red-green-blue
setpixel(scr, x, y, i, i+64, i+128);
					break;
                case 9: // for high amount of iterations
setpixel(scr, x, y, (i >> 11) % 32 << 3, (i >> 5) % 64 << 2, i % 32 << 3);
					break;
				case 10: // same as above, red and blue swapped
setpixel(scr, x, y, i % 32 << 3, (i >> 5) % 64 << 2, (i >> 11) % 32 << 3);
					break;
				default: // grayscale fallback
                	setpixel(scr, x, y, i, i, i);
			}
            zx0+=(xmax-xmin)/scr->width;    // posun na dalszii bod na rzsaadku
        }
		for(uint k = 0; k < step; ++k)	zy0+=(ymax-ymin)/(yend-ystart); // posun
    }
}

/**
 * Saves a PPM of the Screen framebuffer. Adapted from an older homework.
 * Primarily used for debugging and color gradient design.
 * The logo in the manual was generated using this code.
 *
 * @param scr       The Screen structure to dump
 * @param fname     The output file name
 * @param max		Maximum subpixel value to write into the header
 */
void saveppm(Screen* scr, const char* fname, char max){
	FILE *imgf = fopen(fname, "w");
	uchar r = 0; uchar g = 0; uchar b = 0;
	fprintf(imgf, "P6\n%u\n%u\n%hhu\n", scr->width, scr->height, max);
	for(uint i = 0; i < scr->height; ++i){
		for(uint j = 0; j < scr->width; ++j){
			// low quality upsampling from RGB565 to RGB8
			r = (scr->screen[i][j] & 0x001F) << 3; //0000000000011111 > 11111000
			g = (scr->screen[i][j] & 0x07E0) >> 3; //0000011111100000 > 11111100
			b = (scr->screen[i][j] & 0xF800) >> 8; //1111100000000000 > 11111000
			fwrite(&r, sizeof(char), 1, imgf);
            fwrite(&g, sizeof(char), 1, imgf);
            fwrite(&b, sizeof(char), 1, imgf);
		}
	}
	fclose(imgf);
}

/**
 * Takes the center coordinates and a zoom level and updates the delimitaions
 * array with appropriate values. These simple 4 lines of code were
 * procrastinated the most, but given that time to think, the code is tidy.
 * It's also the reason that this code doesn't involve user interaction.
 *
 * @param delim     The delimitations array - xmin, ymin, xmax, ymax
 * @param centerx   The real component of Z in the center of view
 * @param centery   The imaginary component of Z in the center of view
 * @param zoom      The zoom level
 * @param aspect    The aspect ratio
 */
void updatedelim(FLOAT* delim, FLOAT centerx, FLOAT centery,
		FLOAT zoom, FLOAT aspect){
	delim[0] = centerx - aspect/zoom;
    delim[1] = centery - aspect/zoom;
    delim[2] = centerx + 1.0f/zoom;
    delim[3] = centery + 1.0f/zoom;
}

/**
 * Asks user to input most of the values manually. All arguments are addresses
 * as there is no way of returning them all at once and scanf needs addresses.
 *
 * @param x         Real component of Z
 * @param y         Imaginary component of Z
 * @param cx        Real component of C
 * @param cy        Imaginary component of C
 * @param exp       Exponent applied to Z
 */
void promptuser(FLOAT* x, FLOAT* y, FLOAT* cx, FLOAT *cy, char* exp){
	int succ = 1;
	printf("Zx (FLOAT) := "); succ = scanf("%f", x);
    if (!succ) puts("Failed to load a FLOAT into x.");
    printf("Zy (FLOAT) := "); succ = scanf("%f", y);
    if (!succ) puts("Failed to load a FLOAT into y.");
    printf("Cx (FLOAT) := "); succ = scanf("%f", cx);
    if (!succ) puts("Failed to load a FLOAT into cx.");
    printf("Cx (FLOAT) := "); succ = scanf("%f", cy);
    if (!succ) puts("Failed to load a FLOAT into cy.");
	printf("Exponent (char) := "); succ = scanf("%hhd", exp);
    if (!succ) puts("Failed to load a char into exp.");
}

/**
 * Prints a keyboard help in 1 of the 4 variants depending on the given width
 * of the terminal. Not the best looking code, but it helps greatly while
 * the program is running.
 *
 * @param termw     Width of the terminal
 */
void printhelp(ushort termw){
    if (termw > 157) /* 160+ CPL, 4 lines */ puts("\
 [1] Help        [2] Clear       [3] Random C    [4] Random Z    [5] Demo        [6] xxxxxxxxxx  [7] Set values  [8] Reset all   [9] Save PPM    [0] Save TXT\n\
 [Q] Quit        [W] Zoom in     [E] Set zoom    [R] Cy +        [T] Terminal    [Y] Canvas      [U] Zy +        [I] Set iters   [O] Iters +     [P] Redraw\n\
 [A] Reset C     [S] Zoom out    [D] Cx -        [F] Cy -        [G] Cx +        [H] Zx -        [J] Zy -        [K] Zx +        [L] Iters -\n\
 [Z] Palette <   [X] Palette >   [C] Gradient <  [V] Gradient >  [B] Reset Z     [N] Exponent -  [M] Exponent +");
	else if (termw > 117 && termw < 158) /* 120ish CPL, 8 lines */ puts("\
  [1]        [2]        [3]        [4]        [5]        [6]        [7]        [8]        [9]        [0]\n\
 Help       Clear    Random C   Random Z     Demo      xxxxxx   Set values  Reset all  Save PPM   Save TXT\n\
\n\
   [Q]        [W]        [E]            [R]            [T]   [Y]            [U]            [I]        [O]        [P]\n\
              [S]                   [D] [F] [G]                         [H] [J] [K]                   [L]\n\
  Quit       Zoom      Set Zoom       Cx & Cy      Terminal  Canvas       Zx & Zy       Set Iters    Iters     Redraw\n\
\n\
  [A] Reset C    [Z] < Palette > [X]    [C] < Gradient > [V]    [B] Reset Z    [N] - Exponent + [M]");
	else if (termw > 71 && termw < 118) /* 80ish CPL, 10 lines */ puts("\
 [1]  [2]    [3]      [4]   [5]   [6]    [7]     [8]     [9]     [0]\n\
Help Clear RandomC RandomZ Demo xxxxxx SetVal ResetAll SavePPM SaveTXT\n\
\n\
Quit       SetZoom            Term Canvas            SetIters    Redraw\n\
 [Q]   [W]   [E]      [R]      [T]   [Y]      [U]      [I]   [O]   [P]\n\
       [S]         [D][F][G]               [H][J][K]         [L]\n\
      Zoom          Cx & Cy                 Zx & Zy         Iters\n\
\n\
   [A]       [Z]  [X]    [C]  [V]     [B]     [N]  [M]\n\
 ResetC      Palette     Gradient   ResetZ    Exponent");
	else /* 40 CPL, 10 lines */ puts("\
[1][2][3]  [4]  [5] [6] [7] [8] [9] [0]\n\
? Cl RndC RndZ Demo xxx Set Rst PPM TXT\n\
\n\
Quit  SetZoom  TermCnvs   SetIt   Draw\n\
[Q] [W] [E] [R] [T][Y] [U] [I] [O] [P]\n\
    [S]  [D][F][G]  [H][J][K]  [L]\n\
   Zoom   Cx & Cy    Zx & Zy  Iters\n\
\n\
 [A]  [Z][X]  [C][V]  [B]  [N][M]\n\
RstC   Pal     Grad  RstZ   Exp");
}

/**
 * # * *** * J U L I A S C I I * F R A C T A L * *** *
 *
 * Algorithm adapted from [root.cz/clanky/fraktaly-v-pocitacove-grafice-x]
 * (http://www.root.cz/clanky/fraktaly-v-pocitacove-grafice-x/).
 * Formula: Z = Z^n + C, where Z and C are complex numbers and n is an integer.
 *
 * The main function does argument processing, command handling & state holding.
 *
 * @author		getmajan AKA 637man AKA j-61m
 * @version		1.3
 * @date		May 2020
 * @warning  	Makes use of ANSI escape sequences to clear the screen and move
 *              the cursor to top left corner. Windows before 10 support those
 *              only in NTVDM and not in regular Win32 console.
 * @warning     ANSI palette 2 and 3 may not work with your terminal, depending
 *              on the encoding. Should be set to Latin 1 or Western Europe.
 * @copyright	GNU General Public License
 *
 * @param argc	There are 4: buffer input, canvas width and height, and threads.
 *              None is required, but the 2nd isn't accepted without the 3rd
 *              and the 4th isn't accepted without all the previous.
 * @param argv	Implicit values are 1, 480, 320 and 1.
 */
int main(int argc, char** argv){
	/*   ***   I N I T I A L I Z A T I O N   ***   */
	// The order of declarations is optimized
	uint scrx = 480; ///__scrx__: Holds the Screen framebuffer width.
	uint scry = 320; ///__scry__: Holds the Screen framebuffer height.
	FLOAT x = DEFX; ///__x__: Holds the real Z component in center of view.
	FLOAT y = DEFY; ///__y__: Holds the imaginary Z component in center of view.
	FLOAT cx = DEFCX; ///__cy__: Holds the real component of C.
	FLOAT cy = DEFCY; ///__cy__: Holds the imaginary component of C.
    FLOAT delim[4] = {-1.0f,-1.0f,1.0f,1.0f}; ///__delim__: xmin,ymin,xmax,ymax
	FLOAT zoom = DEFZOOM; ///__zoom__: Can zoom up to around 2^18.
    FLOAT aspect = scrx/scry; ///__aspect__: Holds the aspect ratio.
	uint iters = DEFITERS; ///__iters__: Sets the amount of iterations.
    ushort termw = DEFTERMW; ///__termw__: Sets the terminal width.
	ushort termh = DEFTERMH; ///__termh__: Sets the terminal height.
    uchar gradch = DEFGRADCH; ///__gradch__: Chooses the fractal gradient.
	uchar palch = DEFPALCH; ///__palch__: Chooses the fractal ASCII palette.
	char exp = DEFEXP; ///__exp__: Holds the exponent for Z in the formula.
    char redraw = 1; ///__redraw__: Controls redrawing after a command.
    char demo = 0; ///__demo__: Signals to read command from file.
	char drawnow = 1; ///__drawnow__: Signals whether a redraw is necessary.
    char statusbar = 1; ///__statusbar__: Signals whether to print a status bar.
    /**
	 * __command__: Holds the current command to perform. Commands are these:
	 * @code{.unparsed}
	 *  [1]   [2]   [3]      [4]    [5]   [6]    [7]     [8]     [9]     [0]
	 * Help Clear RandomC RandomXY Demo xxxxxx SetVal ResetAll SavePPM SaveTXT
	 *
	 * Quit       SetZoom            Term Canvas           SetIters     Redraw
	 *  [Q]   [W]   [E]      [R]      [T]   [Y]      [U]      [I]   [O]   [P]
	 *        [S]         [D][F][G]               [H][J][K]         [L]
	 *       Zoom          Cx & Cy                  X & Y          Iters
	 *
	 *    [A]       [Z]  [X]    [C]  [V]     [B]     [N]  [M]
	 *  ResetC      Palette     Gradient   ResetZ    Exponent
	 * @endcode
	 */
	char command = 0;
    char scan = 1; ///__scan__: Use scanf(), needs EOLs to not block.
    char udplisten = 0; ///__udplisten__: Listening to datagrams, NYI.
    const uchar PALNUM = 4; ///__PALNUM__: Number of ANSI palettes implemented.
    const uchar GRADNUM = 11; ///__GRADNUM__: Number of gradients implemented.
	int succ = 1; /// __succ__: Used for normally ignored return values.

    uint thrnum = 1; ///__thrnum__: Number of threads.

    if (argc > 1) scan = atoi(argv[1]);
    if (argc > 3){
    	scrx = atoi(argv[2]);
		scry = atoi(argv[3]);
	}
	Screen* scr = newscreen(scrx, scry);

	FILE* demofile = NULL;
	char demofn[64];
	char imagefn[64];
	char textfn[64];

	char* palette[4] = {
	" .:-=+*#%@%#*+=-:.", // http://paulbourke.net/dataformats/asciiart/
	" .`-_':,;^=+/\\\"|)\\<>)iv%xclrs{*}I?!][1taeo7zjLunT#JCwfy325Fp6mqSghVd4Eg\
XPGZbYkOA&8U$@KHDBWNMR0Q", // www.fatvat.co.uk/2009/09/generating-ascii-art.html
	"ANEOOMÉËEAÂWQBAa#NÁ?EÄAHKRŽoXg?eqUŠOÔA€ßpmaâG¶o?é8ÚÜ$ëdUýeÓ?ÖayObYFDnáZPäš\
Çahu§ÝkY®S9žUTe6µOyxÎ3f4o5ôú&aü™2uçw©YL0VÍL±3IIóC@nöoscu‰11‡zJf%¤Itocîrjv1lí=ii\
<>i7†[??×}*{+()\\/»«•¬|!!÷¦——^a„”“~3o2–°­1‹›;:’‘‚’~^¸…·¨´` ",
	// http://ashiskumar.com/create-ascii-gradient/
    " \xb0\xb1\xb2\xdb\xdb\xb2\xb1\xb0 " // block shades
	};
	uchar pallen[4] = {19, 97, 209, 11};

	if (argc > 4) thrnum = atoi(argv[4]);
    pthread_t thrid[thrnum+1]; ///__thrid__: Array of thread IDs.
	Drawargs* drawargs = (Drawargs*)calloc(thrnum,sizeof(Drawargs));
	for(uint i = 0; i < thrnum; ++i){
		drawargs[i].delims = calloc(4,sizeof(FLOAT));
	}

	puts("   J U L I A S C I I   F R A C T A L"); // welcome screen for terminal
    puts("***************************************   coded by 637man 2018/5");

	setdrawargs(drawargs, scr, iters, delim[0], delim[1],
			delim[2], delim[3], cx, cy, exp, gradch, 0, 1);
	drawjulia((void*) drawargs);

	putscreen(scr, termw, termh-12, palette[palch], pallen[palch]);
	printhelp(termw);
	printf("Ready. Enter 1 for help.\n> ");

    while(command != 'q'){
        /*   ***   I N P U T   P R O C E S S I N G   ***   */
		if(scan) command = getchar_unlocked();

		if(demo) succ = fscanf(demofile, "%c", &command);
		if(!succ) puts("Failed to scan a string to demofile.");
	
        switch(command){ // TODO: Make possible for increments to be set
        	case '1': printhelp(termw);drawnow = 0;statusbar = 0;break; // help
            case '2': printf("\x1b[2J\x1b[1;1H"); break; // clear the screen
            case '3': // generate random C
				cx = (FLOAT)random() / RAND_MAX * 4 - 2;
				cy = (FLOAT)random() / RAND_MAX * 4 - 2;
				break;
            case '4': // generate random Z
				x = (FLOAT)random() / RAND_MAX * 4 - 2;
				y = (FLOAT)random() / RAND_MAX * 4 - 2;
				break;
            case '5': // load animation file
                printf("Demo file := "); succ = scanf("%s", demofn);
				if(!succ){
					puts("Failed to scan a string to demofn.");
					break;
				}
				// file can open itself for infinite loop
                if (demofile) fclose(demofile);
              	demofile = fopen(demofn, "r");
				demo = 1; scan = 0;
				break;
			case '6': // reserved (was toggle UDP)
				drawnow = 0; statusbar = 0;
				break;
            case '7': // manually enter values
				promptuser(&x, &y, &cx, &cy, &exp);
				break;
            case '8': // reset everything
            	x = DEFX; y = DEFY; cx = DEFCX; cy = DEFCY;
				zoom = DEFZOOM; iters = DEFITERS;
				demo = 0; udplisten = 0; termw = DEFTERMW; termh = DEFTERMH;
				redraw = 1; palch = DEFPALCH; gradch = DEFGRADCH; exp = DEFEXP;
				break;
            case '9': // save a PPM of the Screen framebuffer
                printf("Image file := "); succ = scanf("%s", imagefn);
                drawnow = 0; statusbar = 0;
				if (!succ) puts("Failed to load a string to imgagefn.");
				else saveppm(scr, imagefn, 255);
				break;
            case '0': // redraw into text file
            	printf("Text file := "); succ = scanf("%s", textfn);
                drawnow = 0; statusbar = 0;
				if (!succ) puts("Failed to load a string to textfn."); else
				savescreen(scr,textfn,scrx,scry,palette[palch],pallen[palch]);
				break;

			case 'q': goto exit; break; // quit
            case 'w': zoom *= 1.25f; break; // zoom in
            case 's': zoom /= 1.25f; break; // zoom out
            case 'e': // set zoom
				printf("Zoom (FLOAT) := "); succ = scanf("%f", &zoom);
    			printf("\x1b[2J\x1b[1;1H");
				if (!succ) puts("Failed to load a FLOAT into zoom.");
				break;
            case 'd': cx -= 0.025f / zoom; break; // decrease real C part
            case 'f': cy -= 0.025f / zoom; break; // decrease imaginary C part
            case 'r': cy += 0.025f / zoom; break; // increase real C part
            case 'g': cx += 0.025f / zoom; break; // increase imaginary C part
            case 't': // set terminal dimension
        		printf("Terminal width (ushort) := ");succ=scanf("%hu",&termw);
                if (!succ) {
					puts("Failed to load an ushort into termw.");
					break;
				}
    			printf("Terminal height (ushort) := ");succ=scanf("%hu",&termh);
    			printf("\x1b[2J\x1b[1;1H");
                if (!succ) puts("Failed to load an ushort into termh.");
				break;
            case 'y': // set canvas dimensions (realloc Screen)
    			printf("Canvas width (ushort) := "); succ = scanf("%u", &scrx);
                if (!succ){
					puts("Failed to load a FLOAT into scrx.");
					break;
				}
    			printf("Canvas height (ushort) := "); succ = scanf("%u", &scry);
                printf("\x1b[2J\x1b[1;1H");
				if (succ){
					freescreen(scr);
					scr = newscreen(scrx, scry);
				}else{
                    if (!succ) puts("Failed to load a FLOAT into scry.");
				}
				break;
            case 'h': x -= 0.05f / zoom; break; // move left
			case 'j': y += 0.05f / zoom; break; // move down (somehow inverted)
            case 'u': y -= 0.05f / zoom; break; // move up (somehow inverted)
            case 'k': x += 0.05f / zoom; break; // move right
            case 'i': // set iteration amount
				printf("Iterations (uint) := "); succ = scanf("%u", &iters);
				printf("\x1b[1;1H");
                if (!succ) puts("Failed to load an uint into iters.");
				break;
            case 'o': iters?iters*=2:++iters; break; // more iterations
            case 'l': iters?iters/=2:--iters; break; // less iterations (or max)
            case 'p': redraw = !redraw; break; // toggle redrawing

            case 'a': cx = DEFCX; cy = DEFCY; break; // reset the C constant
            case 'z': palch = (palch-1) % PALNUM; break; // previous palette
            case 'x': palch = (palch+1) % PALNUM; break; // next palette
            case 'c': gradch = (gradch-1) % GRADNUM; break; // previous gradient
            case 'v': gradch = (gradch+1) % GRADNUM; break; // next gradient
            case 'b': x = DEFX; y = DEFY; break; // reset position
            case 'n': --exp; break; // decrease exponent
            case 'm': ++exp; break; // increase exponent

			case '\n':
            	printf("> ");
				if (scan){
					drawnow = 0;
					statusbar = 0;
				}
				break;
		    case ' ':
            case '\t':
            case '\r': drawnow = 0; statusbar = 0; break; // whitespace
			
			case 255:
		    case -1: // EOF is -1, which is 255 when char is unsigned
			    drawnow = 0; statusbar = 0;
				scan = 1; demo = 0; // demo has ended
				if(demofile) fclose(demofile);
				break;
		    case 0:
		        break;
			default:
				printf("%c (%d) doesn't do anything. Enter 1 for help.\n",
						command, command);
            	drawnow = 0; statusbar = 0;
				//scan = 1; demo = 0; // demo script error (or comment?)
                //if(demofile) fclose(demofile);
		}

        /*   ***   D R A W I N G   ***   */
        updatedelim(delim, x, y, zoom, aspect);
		if(redraw && drawnow){ // Divide et impera.

			for(uint t = 0; t < thrnum; ++t){
            setdrawargs(drawargs+t,scr,iters,delim[0],delim[1],
						delim[2],delim[3],cx,cy,exp,gradch,t,thrnum);
			}
			for(uint t = 0; t < thrnum; ++t){
				pthread_create(thrid + t, NULL, (void*) drawjulia,
						(void*) (drawargs + t));
			}
            for(uint t = 0; t < thrnum; ++t) pthread_join(thrid[t], NULL);

            printf("\x1b[1;1H"); // move cursor to the beginning, Win10+
			putscreen(scr, termw, termh, palette[palch], pallen[palch]);
		}
		drawnow = 1;

		if (statusbar){ // put a status bar at the bottom
        printf("Z:%.3f%+.3fi C:%.3f%+.3fi Zoom:%.2f i:%u n:%hhd Term:%hdx%hd \
Pal:%hhd Grad:%hhd Range:%.3f%+.3fi~%.3f%+.3fi Res:%ux%u \n",
			x, y, cx, cy, zoom, iters, exp, termw, termh, palch, gradch,
			delim[0], delim[1], delim[2], delim[3], scrx, scry);
		}else{
			statusbar = 1;
		}
	}

	exit:
    for(uint i = 0; i < thrnum; ++i) free(drawargs[i].delims);
	free(drawargs);
	freescreen(scr);
    if(demofile) fclose(demofile);
	return 0;
}
