from PIL import ImageDraw
import numpy as np

def ComputeSSD(TODOPatch, TODOMask, textureIm, patchL):
    patch_rows, patch_cols, patch_bands = np.shape(TODOPatch)
    tex_rows, tex_cols, tex_bands = np.shape(textureIm)
    ssd_rows = tex_rows - 2 * patchL
    ssd_cols = tex_cols - 2 * patchL
    assert(ssd_rows > 0)
    assert(ssd_cols > 0)
    # Manipulating arrays for correct types, need floats for image arrays
    TODOPatch = 1.0 * TODOPatch
    textureIm = 1.0 * textureIm
    
    #Converts the mask into a boolean array, entries with 0 become true and all others false
    bool_mask = np.array(TODOMask == 0, dtype='bool')  
    
    SSD = np.zeros((ssd_rows,ssd_cols))
    for r in range(ssd_rows):
        for c in range(ssd_cols):
            # Compute sum square difference between textureIm and TODOPatch
            # for all pixels where TODOMask = 0, and store the result in SSD
            """
            Indexing and using numpy library to efficiently perform the SSD calculation between
            the selected patch in the image and the given TODOPatch
            """
            tex_row = r + 2*patchL
            tex_col = c + 2*patchL
            ssds = np.sum(np.square(TODOPatch[:,:,:] - textureIm[r:tex_row + 1, c:tex_col + 1,:]), axis=2)
            SSD[r, c] = np.sum(ssds, where=bool_mask)
            
        
    return SSD

def CopyPatch(imHole,TODOMask,textureIm,iPatchCenter,jPatchCenter,iMatchCenter,jMatchCenter,patchL, selectPatch):
    patchSize = 2 * patchL + 1
    for i in range(patchSize):
        for j in range(patchSize):
            # Copy the selected patch selectPatch into the image containing
            # the hole imHole for each pixel where TODOMask = 1.
            # The patch is centred on iPatchCenter, jPatchCenter in the image imHole
            # copies pixel-by-pixel into the corresponding indices
            if (TODOMask[i, j] == 1):
                im_row_index = iPatchCenter - patchL + i
                im_col_index = jPatchCenter - patchL + j
                imHole[im_row_index, im_col_index] = selectPatch[i, j]
            
        
    return imHole

def DrawBox(im,x1,y1,x2,y2):
    draw = ImageDraw.Draw(im)
    draw.line((x1,y1,x1,y2),fill="white",width=1)
    draw.line((x1,y1,x2,y1),fill="white",width=1)
    draw.line((x2,y2,x1,y2),fill="white",width=1)
    draw.line((x2,y2,x2,y1),fill="white",width=1)
    del draw
    return im

def Find_Edge(hole_mask):
    [cols, rows] = np.shape(hole_mask)
    edge_mask = np.zeros(np.shape(hole_mask))
    for y in range(rows):
        for x in range(cols):
            if (hole_mask[x,y] == 1):
                if (hole_mask[x-1,y] == 0 or
                        hole_mask[x+1,y] == 0 or
                        hole_mask[x,y-1] == 0 or
                        hole_mask[x,y+1] == 0):
                    edge_mask[x,y] = 1
    return edge_mask
