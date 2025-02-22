import numpy as np
from scipy.optimize import curve_fit

gridSize        = 24
separationStyle = 'hex'
sparsity        = 2  # denotes how close together the squares can be, sparsity of 1 means a chessboard pattern
totalCoords     = 45

_1sqCoords = [147,243,339,435,101,197,293,389,485,151,247,343,439,105,201,
              297,393,489,155,251,347,443,109,205,301,397,493,159,255,351,
              447,113,209,305,401,497,163,259,355,451,117,213,309,405,501]  # not ordered

'''Check projectExperiments.xslx sheet for details'''
patternID = {1	:[101],
             2	:[105],
             3	:[109],
             4	:[113],
             5	:[117],
             6	:[147],
             7	:[151],
             8	:[155],
             9	:[159],
             10	:[163],
             11	:[197],
             12	:[201],
             13	:[205],
             14	:[209],
             15	:[213],
             16	:[243],
             17	:[247],
             18	:[251],
             19	:[255],
             20	:[259],
             21	:[293],
             22	:[297],
             23	:[301],
             24	:[305],
             25	:[309],
             26	:[339],
             27	:[343],
             28	:[347],
             29	:[351],
             30	:[355],
             31	:[389],
             32	:[393],
             33	:[397],
             34	:[401],
             35	:[405],
             36	:[435],
             37	:[439],
             38	:[443],
             39	:[447],
             40	:[451],
             41	:[485],
             42	:[489],
             43	:[493],
             44	:[497],
             45	:[501],
             46	:[209,247,259,301,393],
             47	:[205,251,297,389,447],
             48	:[197,255,347,401,439],
             49	:[201,293,351,355,443],
             50	:[251,305,343,397,451],
             51	:[105,109,113,117,155,159,243,309,343,351,355,405,443,451,485],
             52	:[101,109,117,147,155,197,305,309,339,343,351,401,451,485,497],
             53	:[151,163,197,201,209,213,259,301,339,347,393,401,435,439,489],
             54	:[113,159,205,209,243,251,255,301,347,355,393,405,439,443,447],
             55	:[105,151,163,201,213,247,259,293,297,389,397,435,489,493,501],
             56	:[101,113,147,159,209,243,485],
             57	:[101,113,159,205,209,213,497],
             58	:[147,163,209,243,255,443,485],
             59	:[109,201,247,301,309,355,501],
             60	:[117,151,259,347,351,389,439],
             61	:[163,197,201,247,301,447,501],
             62	:[117,151,259,355,393,405,439],
             63	:[105,117,339,347,351,389,493],
             64	:[147,163,197,255,443,447,501],
             65	:[109,201,309,355,393,405,439],
             66	:[101,205,213,397,401,451,497],
             67	:[155,251,293,297,305,489,493],
             68	:[105,293,297,339,389,489,493],
             69	:[155,251,305,343,435,451,489],
             70	:[305,343,397,401,435,451,497],
             71	:[101,113,147,159,163,197,205,209,213,243,255,443,447,485,501],
             72	:[101,113,147,159,205,209,213,243,255,343,397,401,451,485,497],
             73	:[109,147,163,197,201,209,243,247,255,301,309,443,447,485,501],
             74	:[109,117,151,201,247,259,301,309,347,351,355,389,393,405,439],
             75	:[101,113,155,159,205,213,251,305,343,397,401,435,451,489,497],
             76	:[105,117,155,251,259,293,297,305,339,347,351,389,435,489,493],
             77	:[109,151,163,197,201,247,301,309,355,393,405,439,443,447,501],
             78	:[105,117,151,259,293,297,339,347,351,355,389,393,405,439,493],
             79	:[105,155,251,293,297,305,339,343,397,401,435,451,489,493,497],
             80	:[101,147,197,243,293,339,435],
             81	:[105,151,201,247,297,343,439],
             82	:[109,155,205,251,301,347,443],
             83	:[113,159,209,255,305,351,447],
             84	:[117,163,213,259,309,355,451],
             85	:[101,151,197,293,389,435,485],
             86	:[105,155,201,297,393,439,489],
             87	:[109,159,205,301,397,443,493],
             88	:[113,163,209,305,401,447,497],
             89	:[117,147,213,309,405,451,501],
             90	:[151,247,293,343,389,439,485],
             91	:[155,251,297,347,393,443,489],
             92	:[147,243,309,339,405,435,501],
             93	:[159,255,301,351,397,447,493],
             94	:[163,259,305,355,401,451,497],
             95	:[101,105,147,151,197,201,243,247,293,339,343,389,435,439,485],
             96	:[105,109,155,201,205,247,251,297,301,343,347,393,439,443,489],
             97	:[101,117,147,151,197,213,243,293,309,339,389,405,435,485,501],
             98	:[105,151,155,197,201,247,251,293,297,343,389,393,439,485,489],
             99	:[109,113,159,205,209,255,301,305,347,351,397,401,443,447,493],
             100:[101,117,147,163,213,243,259,309,339,355,405,435,451,497,501],
             101:[109,155,159,205,251,255,297,301,347,351,393,397,443,489,493],
             102:[113,117,163,209,213,259,305,309,355,401,405,447,451,497,501],
             103:[113,159,163,209,255,259,305,351,355,397,401,447,451,493,497],
             104:[105, 163, 255, 347, 401],
             105:[109, 201, 259, 351, 447],
             106:[113, 205, 251, 297, 355],
             107:[155, 209, 301, 393, 443],
             108:[159, 251, 305, 397, 451],
             109:[203, 257, 353, 395, 438],
             110:[200, 230, 298, 306, 342],
             111:[150, 201, 253, 271, 307],
             112:[105, 163, 255, 347, 401],
             999:[101,105,109,113,117,147,151,155,159,163,197,201,205,209,213,
             243,247,251,255,259,293,297,301,305,309,339,343,347,351,355,389,393,
             397,401,405,435,439,443,447,451,485,489,493,497,501]}

def get_patternID(sqSet):
    for k,v in patternID.items():
        if v == sqSet:
            return int(k)

def get_patternIDlist_for_nSq_pattern(patternIDnSq):
    spotlist = []
    spotcoords = patternID[patternIDnSq]
    for spot in spotcoords:
        spotlist.append(get_patternID([spot]))
    return spotlist
            
################# Calibration Details ################

calibration_mapping_file = "\\Lab\\Projects\\EI_Dynamics\\Protocols\\Configurations\\21-12-24_Polygon_Calibration_Map_40x.map"
# 40xWI objective, glass slide, 0.5x camera magnification, polygon numbered grid calibration

'''
Data from 24 Dec 2021 Calibration
    Polygon Frame Fraction		Camera Pixel Number	
    x	    y	                cx	    cy
    0.25	0.25	            285	    341
    0.5	    0.5	                642	    535
    0.75	0.75	            1002	730
    0.125	0.125	            106	    244
    0.625	0.25	            822	    632
    0.375	0.375	            464	    439
    0.875	0.875	            1183	827
'''
map = {
        'x' : [0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875], # fraction polygon frame X
        'y' : [0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875], # fraction polygon frame Y
        'cx': [106,   285,  464,   642, 822,   1002, 1183 ], # camera pixel number X
        'cy': [244,   341,  439,   535, 632,   730,  827  ]  # camera pixel number Y
}

def pixel_scaling(x,m,c):
    y = m*np.array(x) + c
    return y

    
poptx,_  = curve_fit(pixel_scaling,map['x'],map['cx'])
popty,_  = curve_fit(pixel_scaling,map['y'],map['cy'])

x0,y0    = [pixel_scaling(0,*poptx), pixel_scaling(0,*popty)]
x1,y1    = [pixel_scaling(1,*poptx), pixel_scaling(1,*popty)]

polygon_frame_properties =  {   
                            'top_left'      :[x0,y0],
                            'top_right'     :[x1,y0],
                            'bottom_right'  :[x1,y1],
                            'bottom_left'   :[x0,y1],
                            'width um'      : round(x1-x0),
                            'height um'     : round(y1-y0),
                            'aspect_ratio'  :(x1-x0) / (y1-y0),
                            'scaling'       :[round((x1-x0)/gridSize), round((y1-y0)/gridSize)],
                            'offsetx'       : round(x0),
                            'offsety'       : round(y0),
                            }

polygon_protocol_patterns_per_sweep_LUT = {
                                        '2_210303_24hex_1sq_ExtFreq_3repeats_135frames' : 1,
                                        '3_210303_24hex_5-15sq_ExtFreq_3repeats_30frames'   : 1,
                                        '3_210331_24hex_5-15sq_ExtFreq_3repeats_24frames'   : 1,
                                        '3_210428_24hex_7-15sq_LTM_rand_ExtFreq_24frames'   : 1,
                                        '3_210428_24hex_7-15sq_LTM_rand_ExtFreq_3repeats_24frames'  : 1,
                                        '3_210428_24hex_7-15sq_LTM_Seq_ExtFreq_24frames'    : 1,
                                        '3_210428_24hex_7-15sq_LTM_Seq_ExtFreq_3repeats_24frames'   : 1,
                                        '3_220117_24hex_7-15sq_LTM_rand_ExtFreq_3repeats_72frames'  : 1,
                                        '5-210723_24hex_15sq_Convergence_IntFreqExtFrame_2ms_1repeat_18frames_3patternperSweep' : 3,
                                        '5-210723_24hex_7-15sq_Convergence_IntFreqExtFrame_2ms_1repeat_36frames_3patternperSweep'  : 3,
                                        '5-210723_24hex_7sq_Convergence_IntFreqExtFrame_2ms_1repeat_18frames_3patternperSweep'  : 3,
                                        '6_221107_24hex_15sq_Convergence_ExtFreq_1repeat_24frames'  : 8,
                                        '6_221107_24hex_7sq_Convergence_ExtFreq_1repeat_24frames'   : 8,
                                        '7_221108_24hex_15sq_Convergence+PulseTrain_ExtFreq_1repeat_8sweeps'    : 20,
                                        '7_221108_24hex_5sq_Convergence+PulseTrain_ExtFreq_1repeat_8sweeps' : 20,
                                        '9_230414_24hex_3sq_Surprise_ExtFreq_10repeats_33frames' :33,
                                        '9_230203_24hex_5sq_Surprise_ExtFreq_10repeats_33frames' :33,
                                        '9_230203_24hex_15sq_Surprise_ExtFreq_10repeats_33frames' :33,
                                        'all_24hex_grid_squares' : 1}


def polygon_protocol_sweep_division(coordfile):
    print(f'coordfile: {coordfile}')
    pattern_per_sweep = polygon_protocol_patterns_per_sweep_LUT[coordfile]
    return pattern_per_sweep

blackfly_camera_pixel_size_um = {'4x': 2.33, '40x': 0.233} # um/pixel
blackfly_camera_resolution_px_per_um = {'4x': 0.429, '40x': 4.29} # pixel/um

polygon_frame_width_px = {'40x': polygon_frame_properties['width um'] * blackfly_camera_resolution_px_per_um['40x'],
                            '4x': polygon_frame_properties['width um'] * blackfly_camera_resolution_px_per_um['4x']
                        }

if __name__ == "__main__":
    for prop,prop_val in polygon_frame_properties.items():
        print(f'{prop} : {prop_val}')