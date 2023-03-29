import sys
from helpers import set_actives_dep_weights, set_target_height, percent_to_height
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
            
        self.degree        = degree
        self.rows          = rows
        self.cols          = cols
        self.numUmbrellas  = self._get_num_cells()
        self.center_pos    = self._get_center_position()
        self.num_links     = self._get_num_links()
        self.min_height    = min_height
        self.height_scales = height_fct(self.numUmbrellas) if height_fct is not None else [1]*self.numUmbrellas
    
    def generate_mesh(self, json_filename, verbose=True):
        genUmbrellaWithHeights(self.degree,
                               self.rows,
                               self.cols,
                               self.height_scales,
                               self.min_height,
                               json_filename=json_filename)
        input_path = f'../UmbrellaGen/{json_filename}.json.gz'
        
        # consant plate height: should check input_data['bbox_diag']
        io, self.input_data, target_mesh, self.curr_um, plate_thickness_scaled, target_height_multiplier\
            = parse_input(input_path,
                          isHex = (self.degree==6),
                          use_target_surface = False)
        
        self.init_heights    = self.curr_um.umbrellaHeights
        self.plate_thickness = self.input_data['thickness']
        
        if verbose:
            print(f"PLATE CHARACTERISTIQUES:\n\
\tplate thickness   : {self.plate_thickness:.6f}, {plate_thickness_scaled}\n\
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
            print('success : ', success, '\n')
            print('energies:')
            print(*allEnergies(self.curr_um).items(), sep='\n')
        if show_plots:
            eqays.plot()

    def border_cells(self):
        if self.degree==3:
            step = 2*self.cols
            if self.rows==1:
                return [0, step-1]
            
        if self.degree==4:
            step = self.cols
        bot_left  = 0
        bot_right = step-1
        top_left  = step*(self.rows-1)
        top_right = step*self.rows
        return  list(range(bot_right)) + \
                list(range(bot_right, top_right-1, step)) + \
                list(range(step, top_left, step)) + \
                list(range(top_left,top_right))

    def center_cells(self):
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
            nb_links += self.cols*(2*self.rows-1)
        elif self.degree==4:
            nb_links += self.cols*(  self.rows-1)
        return nb_links
    
     def _get_center_position(self):
        center_position = np.zeros([self.numUmbrellas, 3])
        for i in range(self.numUmbrellas):
            top_idx = self.curr_um.getUmbrellaCenterJi(i, 0)
            center_position[i] = self.curr_um.joint(top_idx).position
        return center_position
        