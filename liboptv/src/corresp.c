/*  finding correspondences between 2D images of particles, for the purpose
    of composing 3D positions from them.
    
    Algorithm based on the one from 3D-PR\TV, but generalized to any number
    of cameras.
*/

#include "tracking_frame_buf.h"
#include "epi.h"

/* qs_coord2d_x() is a helper function doing inplace quicksort of a 2D
   coordinates array. The sort key is the X coordinate.
   
   Arguments:
   coord_2d *crd - the array of coordinates to sort.
   int left, right - indices of the edges of the currently-sorted part of the
      array (first call should use 0 and length-1).
*/
void qsort_coord2d_x(coord_2d *crd, int left, int right) {
    int over, under;
    double median;
    coord_2d temp;
    
    /* initialize switch-pair search */
    over = left; 
    under = right;
    median = crd[(left + right)/2].x;
    
    do {
        /* Left point over median and right point under median - switch. */
        while (crd[over].x < meadian  &&  over < right)
            over++;
        while (median < crd[under].x  &&  under > left)
            under--;

        if (over <= under) {
            temp = crd[over];
            crd[over] = crd[under];
            crd[under] = temp;
            
            over++;
            underj--;
        }
    }
    while (over <= under);
}


#define pair_ix(p1, p2, num_cams) ((p1)*(num_cams) + (p2))

/* correspondences() generates a list of correspondence groups. Each
   correspondence group is the set of target numbers identifying each one
   target in one image from the set, such that each target in the cgroup is an
   image of the same 3D particle.

   Arguments:
   frame *frm - holds all needed image data.
   calibration **calib - holds camera calibration information needed for
       getting real metric target coordinates corrected from those identified 
       in image. Array of num_cams pointers to calibration objects.
   volume_par *vpar - search volume parameters, contains the search volume for
       epipolar lines and the measure of correspondence needed for a candidate
       target.
   
   Returns:
   n_tupel **corres_lists - a num_cams list of corres-lists, each is a list of 
        4-item arrays containing the index of target for each particle in the
        list in each of up to 4 images, in order. Each list is terminated by an
        n_tupel with negative .corr and garbage p[]. The top-level array is by
        number of particles in correspondence, so that the first element is
        always NULL and is there for convenience only.
        
        Although internally the code uses x-sorting of targets, returned target 
        indices (pnr) are into the *unsorted* target array, as passed in by the
        user.
*/
n_tupel **correspondences(frame *frm, Calibration **calib, volume_par *vpar,
    control_par *cpar)
{
    int cam, part, pair, num_cams, num_pairs;
    int part_img, epi_img;
    double img_x, img_y;
    double epi_start[2], epi_end[2]; /* each for x,y coordinates of a point */
    target *targ; /* shortcut to working target. */
    int num_cands, cix;
    coord_2d **corrected;
    correspond *img_corr;
    correspond **list; /* per-camera-pair list of correspondence candidates
                          lists (one per master target) */
    int *cam_set, *cam_subset, *subset_count, *subset_iter;
    int subset, subset_size, pow_set_size;
    
    num_cams = frm->num_cams; /* save dereferences */
    num_pairs = num_cams * (num_cams - 1) / 2;
    
    /* We work on distortion-corrected image coordinates of particles.
       This loop does the correction. It also recycles the iteration on
       frm->num_cams to allocate some arrays needed later and do some related
       preparation. */
    corrected = (coord_2d **) malloc(num_cams * sizeof(coord_2d *));
    for (cam = 0; cam < num_cams; cam++) {
        corrected[cam] = (coord_2d *) malloc(
            frm->num_targets[cam] * sizeof(coord_2d));
        
        for (part = 0; part < frm->num_targets[cam]; part++) {
            pixel_to_metric(&img_x, &img_y,
                frm->targets[cam][part].x, frm->targets[cam][part].y, cpar);
            
            img_x -= calib[cam]->int_par.xh;
            img_y -= calib[cam]->int_par.yh;

            correct_brown_affin (img_x, img_y, calib[cam]->added_par,
               &corrected[cam][part].x, &corrected[cam][part].y);
            
            /* Ensure that after sorting, there is a mapping to the old 
               indices from the sorted array */ 
            corrected[cam][part].pnr = frm->targets[cam][part].pnr;
        }
        qsort_coord2d_x(corrected, 0, frm->num_targets[cam]);
    }
    
    /* Build pairwise correspondence lists of particles in all images to each-
       other, using a list of candidates on the epipolar line from particle in
       one image to a paired image. A link exists if a cell's .corr attribute 
       != 0.
    */
    list = (list **) malloc(num_pairs * sizeof(list *));
    for (part_img = 0; part_img < num_cams - 1; part_img++) {
        for (epi_img = part_img + 1; epi_img < num_cams; epi_img++) {
            /* Note that we always take the particle from the first image
               and candidates from the second. The reverse may give slightly
               different results. However, it is costly to try both and not
               a big improvement. */
            
            img_corr = (list *) malloc(frm->num_targets[part_img] * 
                sizeof(correspond));
            
            for (part = 0; part < frm->num_targets[part_img]; part++) {
                /* Find epipolar line on corrected image */
                epi_mm(corrected[part_img][part].x, 
                       corrected[part_img][part].y,
                       calib[part_img], cpar->mm, vpar, 
                       &(epi_start[0]), &(epi_start[1]), 
                       &(epi_end[0]), &(epi_end[1])
                );
                
                /* Find candidates close to epipolar line */
                targ = &(frm->targets[part_img][corrected[part_img][part].pnr]);
                num_cands = find_candidate(corrected[epi_img], 
                    frm->targets[epi_img], frm->num_targets[epi_img], 
                    epi_start[0], epi_start[1], epi_end[0], epi_end[1], 
                    targ->n, targ->nx, targ->ny, targ->sumg, 
                    cand, vpar, cpar, calib[epi_img]
                );
                
                /* Copy candidate information to the pairwise list. */
                if (num_cands > maxcand) num_cands = maxcand;
                for (cix = 0; cix < num_cands; cix++) {
                    img_corr[part].p2[cix] = cand[cix].pnr;
                    img_corr[part].corr[cix] = cand[cix].corr;
                    img_corr[part].dist[cix] = cand[cix].tol;
                }
                img_corr[part].n = num_cands;
            }
            pair = pair_ix(part_img, epi_img, num_cams);
            list[pair] = img_corr;
        }
    }
    
    /* Match candidates between sets of (n > 1) images. Iteration is over the
       subset of the power-set of 0 .. num_cams - 1, i.e. all combinations of
       2, ..., num_cams members of the camera indices set (unordered
       combinations).
    */
    cam_subset = (int *) malloc(num_cams * sizeof(int));
    subset_iter = (int *) malloc((num_cams - 1) * sizeof(int));
    subset_cand_iter = (int *) malloc((num_cams - 1) * sizeof(int));
    
    pow_set_size = (1 << num_cams) - 1;
    for (subset = 1; subset <= pow_set_size; subset++) {
        /* Prepare the subset for iterating over its particle groups. */
        
        /* popcount, done slowly bcz other stuff also happens */
        subset_size = 0; 
        for (cam = 0; cam < num_cams; cam++) {
            if (subset & cam_set[cam]) {
                cam_subset[subset_size] = cam;
                subset_iter[subset_size] = 0;
                subset_cand_iter[subset_size] = 0;
                subset_size++;
            }
        }
        subset_iter[subset_size - 1] = -1; /* no-candidates case */
        
        /* Go over particles group belonging to this camera subset. */
        carry = 0;
        while (carry == 0) {
            /* Advance to next candidate group*/
            carry = 1;
            
            is_clique = 0;
            for (part_img = 0; part_img < subset_size - 1; part_img++) {
                for (epi_img = part_img + 1; epi_img < subset_size; epi_img++) {
                    pair = pair_ix(cam_subset[part_img], 
                        cam_subset[epi_img], num_cams);
                }
            }
        }
        /***********/
    }
}