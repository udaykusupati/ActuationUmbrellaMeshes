import sys
import numpy as np
from helpers import set_actives_dep_weights, set_target_height, percent_to_height, get_center_position
sys.path.append('..')
from configuration import parse_input, deploy_umbrella_pin_rigid_motion
from visualization_helper import get_color_field
from pipeline_helper import allEnergies
sys.path.append('../UmbrellaGen')
from grid_gen import genUmbrellaWithHeights
class UmbrellaGrid:
    def __init__(self, degree=3, rows=2, cols=2, height_fct=None, min_height=64):
        assert degree==3 or degree==4 or degree==6, f'degree is {degree}, but is should either be 3, 4 or 6'
        if degree==4 and (rows<2 or cols<2):
            raise ValueError('rows and/or cols sould be greater than two.')
            
        self.degree          = degree
        self.rows            = rows
        self.cols            = cols
        self.numUmbrellas    = self._get_num_cells()
        self.num_links       = self._get_num_links()
        self.min_height      = min_height
        self.height_scales   = height_fct(self.numUmbrellas) if height_fct is not None else [1]*self.numUmbrellas
    
    def generate_mesh(self, json_filename, verbose=True):
        genUmbrellaWithHeights(self.degree,
                               self.rows,
                               self.cols,
                               self.height_scales,
                               self.min_height,
                               json_filename=json_filename)
        input_path = f'../UmbrellaGen/{json_filename}.json.gz'
        
        _, self.input_data, _, self.curr_um, _, _ = parse_input(input_path,
                                                                isHex = (self.degree==6),
                                                                use_target_surface = False)
        self.init_heights    = self.curr_um.umbrellaHeights
        self.plate_thickness = self.input_data['thickness']
        self.init_center_pos = get_center_position(self.curr_um)
        
        if verbose:
            print(f"PLATE CHARACTERISTIQUES:\n\
\tplate thickness   : {self.plate_thickness:.6f}\n\
\tplate edge length : {self.input_data['plate_edge_length']:.6f}")
    
    def deploy(self, active_cells, target_percents, view=None, rod_colors=None, uidBased=False, verbose=True, show_plots=False):
        dep_weights              = set_actives_dep_weights(self.numUmbrellas, active_cells)
        target_heights           = percent_to_height(self.init_heights, self.plate_thickness, active_cells, target_percents)
        target_height_multiplier = set_target_height(self.numUmbrellas, active_cells, target_heights)
        if rod_colors is None:
            rod_colors = get_color_field(self.curr_um, self.input_data, uidBased)
        success, eqays = deploy_umbrella_pin_rigid_motion(self.curr_um,
                                                          self.plate_thickness,
                                                          target_height_multiplier,
                                                          view,
                                                          rod_colors,
                                                          analysis = True,
                                                          dep_weights=dep_weights)
        if verbose:
            print('success:', success, '\n')
            ener = 'energies:\n'
            for e in allEnergies(self.curr_um).items():
                ener += f'  {e[1]: 10f}: {e[0]}\n'
            print(ener)
        if show_plots:
            eqays.plot()

    def border(self):
        if self.degree==3:
            step = 2*self.cols
            if self.rows==1:
                return [0, step-1]
            if self.cols==1:
                return [0, 2*self.rows-1]
            
        if self.degree==4: # can't have unitary cols/rows
            step = self.cols
        bot_left  = 0
        bot_right = step-1
        top_left  = step*(self.rows-1)
        top_right = step*self.rows
        return  list(range(bot_right)) + \
                list(range(bot_right, top_right-1, step)) + \
                list(range(step, top_left, step)) + \
                list(range(top_left,top_right))

    def center(self):
        if self.degree==3:
            if self.rows%2==0:
                a = (self.rows-1)*self.cols-1
                b = (self.rows+1)*self.cols
                return [a, a+1, b-1, b]
            else:
                a = self.rows*self.cols
                return [a-1, a]
        if self.degree==4:
            if self.rows%2==0:
                a = int((self.rows/2-0.5)*self.cols)
                b = int((self.rows/2+0.5)*self.cols)
                if self.cols%2==0:
                    return [a-1, a, b-1, b]
                else:
                    return [a, b]
            else:
                a = int((self.rows/2)*self.cols)
                if self.cols%2==0:
                    return [a-1, a]
                else:
                    return [a]
    def vline(self, i=0):
        e_msg = f'vertical line index ({i}) should not be greater than'+ ' {}'
        if self.degree==3: return self._vline(i, e_msg, 2)
        if self.degree==4: return self._vline(i, e_msg, 1)
    

    def hline(self, i=0):
        e_msg = f'horizontal line index ({i}) should not be greater than'+' {}'
        if self.degree==3: return self._hline(i, e_msg, 2)
        if self.degree==4: return self._hline(i, e_msg, 1)
    
    def cross(self, length=0):
        # raise Error is incorect value
        cross = [] # should we get `larger` as parameter ?
        center = self.center_cells()
        for l in range(length):
            for i,c in enumerate(center):
                if self.degree==4:
                    cross.extend([c-l, c+l]) # horizontal
                    cross.extend([c-l*self.cols, c+l*self.cols]) # vertical
                if self.degree==3:
                    a = c-l*(2*self.cols)
                    b = c+l*(2*self.cols)
                    cross.extend([c-l-1, c+l+1]) # horizontal
                    cross.extend([a, b]) # vertical
        return list(np.unique(np.array(cross)))
    
    # =====
    # Helpers
    # =====
    def _get_num_cells(self):
        if self.degree == 3:
            return self.rows*(self.cols*2)
        else:
            return self.rows*self.cols
    
    def _get_num_links(self):
        nb_links = self.rows*(self.cols-1)
        if self.degree==3:
            return nb_links + self.cols*(2*self.rows-1)
        elif self.degree==4:
            return nb_links + self.cols*(  self.rows-1)
        else: raise ValueError('undefined function `_get_num_links` for this degree value')
        
    # V-line
    def _vline(self, i, msg, step=1):
        if i > step*self.cols-1: raise ValueError(e_msg.format(step*self.cols-1))
        line = []
        for j in range(self.rows):
            line+=[i+j*step*self.cols]
        return line
    
    # H-line
    def _hline(self, i, msg, step=1):
        if i > self.rows-1: raise ValueError(msg.format(self.rows-1))
        line = []
        for j in range(step*self.cols):
            line+=[(step*self.cols)*i+j]
        return line
        
