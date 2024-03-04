// Author: APD team, except where source was noted

#include "helpers.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include "pthread_barrier_mac.h"


#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }

// All the threads will share this info
typedef struct General_info {
    uint8_t rescale; // set to 1 if rescale is needed and 0 if not
    int P; // number of threads runnnig in parallel
    ppm_image *image; // initial image
    ppm_image *rescaled_image; // this is only used if the initial image needs rescaling
    int step_x; 
    int step_y;
    unsigned char sigma; // threshold for the Marching Squares Algorithm
    unsigned char **grid; 
    ppm_image **contour_map;
    pthread_barrier_t *barrier;
} general_info;

typedef struct {
    int thread_id;
    general_info *general_info;
} thread_arg;

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "../checker/contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);

    free(image->data);
    free(image);
}

void *thread_function(void *thread_function_argument) {

    thread_arg *arg = (thread_arg*)thread_function_argument;
    general_info *info = (general_info *)arg->general_info;
    ppm_image *rescaled_image = info->rescaled_image;

    // Each separate thread will operate and modify the image only between start & end 
    // p & q are the number of partitions which will result after dividing 
    // the width and height of the image to the preset STEP
    int p, q, start, end;

    if (info->rescale == 1) {
        p = rescaled_image->x / info->step_x;
        q = rescaled_image->y / info->step_y;
    } else {
        p = info->image->x / info->step_x;
        q = info->image->y / info->step_y;
    }

    // 1. (Optional) Check if image needs rescaling
    if (info->rescale) {
        uint8_t sample[3];

        start = arg->thread_id * (double)rescaled_image->x / info->P;
        end = fmin((arg->thread_id + 1) * (double)rescaled_image->x / info->P, rescaled_image->x);

        // Use bicubic interpolation for scaling
        for (int i = start; i < end; i++) {
            for (int j = 0; j < rescaled_image->y; j++) {
                float u = (float)i / (float)(rescaled_image->x - 1);
                float v = (float)j / (float)(rescaled_image->y - 1);
                sample_bicubic(info->image, u, v, sample);

                rescaled_image->data[i * rescaled_image->y + j].red = sample[0];
                rescaled_image->data[i * rescaled_image->y + j].green = sample[1];
                rescaled_image->data[i * rescaled_image->y + j].blue = sample[2];
            }
        }
        // To avoid race conditions, all the threads will wait for each other after every step of the algorithm
        pthread_barrier_wait(info->barrier);
        arg->general_info->image = rescaled_image;
        
    }

    // 2. Sample the grid 
    start = arg->thread_id * (double)p / info->P;
    end = fmin((arg->thread_id + 1) * (double)p / info->P, p);

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = info->image->data[i * info->step_x * info->image->y + j * info->step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > info->sigma) {
                info->grid[i][j] = 0;
            } else {
                info->grid[i][j] = 1;
            }
        }
    }

    // Last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    for (int i = start; i < end; i++) {
        ppm_pixel curr_pixel = info->image->data[i * info->step_x * info->image->y + info->image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > info->sigma) {
            info->grid[i][q] = 0;
        } else {
            info->grid[i][q] = 1;
        }
    }

    // Last thread will manipulate the last column of pixels
    // The following code section will be executed
    // only by the thread with the thread_id = P-1 
    if (arg->thread_id == info->P - 1) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = info->image->data[(info->image->x - 1) * info->image->y + j * info->step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > info->sigma) {
                info->grid[p][j] = 0;
            } else {
                info->grid[p][j] = 1;
            }
        }
    }

    // 3. March the squares

    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            unsigned char k = 8 * info->grid[i][j] + 4 * arg->general_info->grid[i][j + 1] + 2 * info->grid[i + 1][j + 1] + 1 * info->grid[i + 1][j];
            update_image(info->image, info->contour_map[k], i * info->step_x, j * info->step_y);
        }
    }
    pthread_barrier_wait(info->barrier);
    free(arg);
    return NULL;
}

int main(int argc, char *argv[]) {
    // Check correct usage of the program
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    uint8_t rescale = 0; // will be set to 1 if the image needs rescaling
    int P = atoi(argv[3]); // number of parallel running threads
    ppm_image *rescaled_image = NULL;
    int arguments[P];
    for (int i = 0; i < P; i++) {
        arguments[i] = i;
    }

    ppm_image *image = read_ppm(argv[1]);
    int step_x = STEP;
    int step_y = STEP;

    // 0. Initialize contour map
    ppm_image **contour_map = init_contour_map();

    // 1. Check if the image needs rescaling. We only rescale downwards
    if (image->x > RESCALE_X || image->y > RESCALE_Y) {
        rescale = 1; // rescale is needed
        rescaled_image = (ppm_image *)malloc(sizeof(ppm_image));
        if (!rescaled_image) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
        rescaled_image->x = RESCALE_X;
        rescaled_image->y = RESCALE_Y;

        rescaled_image->data = (ppm_pixel*)malloc(rescaled_image->x * rescaled_image->y * sizeof(ppm_pixel));
        if (!rescaled_image) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    // 2. Create thread_argument struct
    int i;
	pthread_t tid[P];
    pthread_barrier_t *barrier = (pthread_barrier_t *)malloc(sizeof(pthread_barrier_t));
    if (!barrier) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }
	pthread_barrier_init(barrier, NULL, P);
    
    int p = rescale == 1 ? RESCALE_X / step_x : image->x / step_x;
    int q = rescale == 1 ? RESCALE_Y / step_y : image->y / step_y;

    unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    for (int i = 0; i <= p; i++) {
        grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    // All the threads will share this informations
    struct General_info *general_info = (struct General_info *)malloc(sizeof(struct General_info));
    if (!general_info) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    // This information will be shared by all the threads
    general_info->rescale = rescale;
    general_info->P = P;
    general_info->image = image;
    general_info->step_x = step_x;
    general_info->step_y = step_y;
    general_info->sigma = SIGMA;
    general_info->grid = grid;
    general_info->contour_map = contour_map;
    general_info->barrier = barrier;
    general_info->rescaled_image = rescaled_image;

	// 3. Start threads
	for (i = 0; i < P; i++) {
        thread_arg *arg = (thread_arg*)malloc(sizeof(thread_arg));
        if (!arg) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
        arg->general_info = general_info;
        arg->thread_id = arguments[i]; 
		pthread_create(&tid[i], NULL, thread_function, arg);
	}

	for (i = 0; i < P; i++) {
		pthread_join(tid[i], NULL);
	}

    // 4. Write output
    write_ppm(general_info->image, argv[2]);

    if (rescale) {
         free(rescaled_image->data);
        free(rescaled_image);
    }
    free(general_info);

    pthread_barrier_destroy(barrier);

    return 0;
}
