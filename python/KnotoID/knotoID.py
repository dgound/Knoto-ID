#%%
import numpy as np
import pandas as pd
import random
import time
import json
import knotoID_cpp as kn
from importlib.resources import files, as_file
import matplotlib.pyplot as plt
import os
import re
import contextlib
#%%
class KnotoID:
    def __init__(self, 
                 filename:str=None, 
                 seed:int=None, 
                 nb_projections:int=None, 
                 projections_list:str=None, 
                 precomputed_projections:np.ndarray=None,
                 closure_method:str='open',
                 planar:bool=False,
                 arrow_polynomial:bool=False,
                 cyclic:bool=False,
                 debug:bool=False):

        """
        Initialize Knoto-ID with given parameters and validate them.
        
        :param projection: Tuple (dx, dy, dz) for fixed projection direction, or None for random.
        :param seed: Seed for random number generator.
        :param planar: Boolean indicating if the diagram is on a plane.
        :param arrow_polynomial: Evaluate arrow polynomial for knotoids if True.
        :param debug: Enable debug mode.
        :param simplify_diagram: Simplify the knot(oid) diagram if True.
        :param max_moves_III: Maximum number of random Reidemeister move III iterations.
        :param max_unsuccessful_moves_III: Maximum number of unsuccessful iterations for move III.
        :param timeout: Time limit for polynomial evaluation in seconds.
        """      


        self.filename = filename
        self.nb_projections = nb_projections
        self.projections_list = projections_list
        self.precomputed_projections = precomputed_projections
        self.flag_planar = planar
        self.flag_arrow_polynomial = arrow_polynomial
        self.closure_method = closure_method
        self.flag_debug = debug
        self.flag_cyclic = False
        self.names_db = None
        self.curve_analysis = None
        self.knotted_core_analysis = None
        self.all_chains_analysis = None
        self.kc_search_path_analysis = None
    
        if seed is not None:
            kn.set_seed(seed)
        else:
            kn.set_seed(int(time.time()))

        if self.projections_list is not None: 
            if self.projections_list == 'internal':
                self.projections, self.weights = self.load_resources('projections')
            elif self.projections_list == 'precomputed':
                self.projections, self.weights = self.precomputed_projections[:,0:3], self.precomputed_projections[:,3]
        else:
            self.projections = [(random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)) for _ in range(self.nb_projections)]
            self.weights = [1.0 / self.nb_projections] * self.nb_projections

        self.polygons = []

        self. validate_class_parameters(projections_list, nb_projections, precomputed_projections,seed,debug)

    def validate_class_parameters(self,projections_list, nb_projections, precomputed_projections,seed,debug):
       
        if projections_list == 'precomputed' and precomputed_projections is None:
            raise ValueError("Precomputed projections must be provided with the precomputed_projections parameter.")

        if projections_list == 'precomputed' and (not isinstance(precomputed_projections, np.ndarray) or precomputed_projections.ndim != 2 or precomputed_projections.shape[1] != 4):
            raise ValueError("Projections list must be a NumPy array of shape (n, 4). The first three columns are the projection directions (dx, dy, dz) and the fourth column is the weight (w).")
        
        if projections_list is None and (nb_projections is None or nb_projections < 1):
            raise ValueError("Number of projections must be a positive integer.")
        
        if seed is not None and not isinstance(seed, (int, float)):
            raise ValueError("Seed must be an integer or None.")
        
        if not isinstance(debug, bool):
            raise ValueError("Debug must be a boolean value.")

        if self.flag_cyclic or self.closure_method =='direct' or self.closure_method =='rays':
            self.names_db = self.load_resources('knot_names')

        elif self.closure_method =='open':
            if self.flag_planar:
                if self.flag_arrow_polynomial:
                    self.names_db = self.load_resources('knotoid_names_planar_arrow')
                else:
                    self.names_db = self.load_resources('knotoid_names_planar')
            else:
                if self.flag_arrow_polynomial:
                    self.names_db = self.load_resources('knotoid_names_sphere_arrow')
                else:
                    self.names_db = self.load_resources('knotoid_names_sphere')
        # Validate closure_method
        if self.closure_method not in ["open", "direct", "rays"]:
            raise ValueError("closure_method must be 'open', 'direct', or 'rays'.")
        
        if self.closure_method=='direct' or self.closure_method=='rays':
            self.flag_cyclic = True
        # if self.flag_cyclic==False and self.closure_method=='direct':
        #     raise ValueError("Closure_method must be 'open' for non-cyclic curves (flag_cyclic==False).")

        # if self.flag_cyclic==False and self.closure_method=='rays':
        #     raise ValueError("Closure_method must be 'open' for non-cyclic curves (flag_cyclic==False).")

        if self.closure_method=='direct':
            # print('Cyclic curve: Forcing nb-projections=1.')
            self.nb_projections = 1
            self.projections = [(random.uniform(-1, 1), random.uniform(-1, 1), random.uniform(-1, 1)) for _ in range(self.nb_projections)]
            self.weights = [1.0 / self.nb_projections] * self.nb_projections

            if self.projections_list is not None:
                # print('Cyclic curve: projections_list is ignored.')
                self.projections_list = None

        #Validate Boolean flags
        if not isinstance(self.flag_planar, bool):
            raise ValueError("planar must be a boolean value.")
        
        if not isinstance(self.flag_arrow_polynomial, bool):
            raise ValueError("arrow_polynomial must be a boolean value.")
        

    def load_resources(self, data: str):
        resources = {
            'projections': 'data/projections_list_100.txt',
            'knot_names': 'data/knot_names.txt',
            'knotoid_names_planar_arrow': 'data/knotoid_names_planar_arrow.txt',
            'knotoid_names_planar': 'data/knotoid_names_planar.txt',
            'knotoid_names_sphere_arrow': 'data/knotoid_names_sphere_arrow.txt',
            'knotoid_names_sphere': 'data/knotoid_names_sphere.txt'
        }


        resource = files('KnotoID')
        for part in resources[data].split('/'):
            resource = resource / part

        with as_file(resource) as file_path:
            if data == 'projections':
                return self.load_projections(str(file_path))

            return pd.read_csv(file_path, names=['Name', 'Polynomial'], sep="\t")


    def load_projections(self,filename):
        projections = []
        weights = []
        total_weight = 0

        with open(filename, 'r') as file:
            for line in file:
                if '#' in line:
                    line = line[:line.find('#')]
                parts = line.strip().split()
                if len(parts) >= 3:
                    x, y, z = map(float, parts[:3])
                    weight = float(parts[3]) if len(parts) == 4 else 1
                    projections.append((x, y, z))
                    weights.append(weight)
                    total_weight += weight

        # Normalize weights
        weights = [w / total_weight for w in weights]
        return np.array(projections), np.array(weights)

    @contextlib.contextmanager
    def _suppress_cpp_output(self, debug=None):
        """Silence the C++ stdout/stderr that the bindings route into Python's
        streams (via the pybind ostream_redirect guards), unless debug is on."""
        active_debug = self.flag_debug if debug is None else debug
        if active_debug:
            yield
            return
        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
                yield

    def update_histogram(self, histogram, polynomial, weight):
        """
        Update the histogram with the given polynomial and its weight.
        """
        if polynomial not in histogram:
            histogram[polynomial] = 0
        histogram[polynomial] += weight
    


    def create_histogram_dataframe(self, histogram,knotted_core):
        """
        Convert the histogram into a Pandas DataFrame, include names from names_db, 
        and sort it in descending order of frequency.

        Args:
        histogram (dict): The histogram to be converted.
        names_db (pd.DataFrame, optional): A DataFrame of names corresponding to polynomials.
        """
        histogram_data = []
        for polynomial_string, info in histogram.items():
            polynomial_string = polynomial_string.replace('<Polynomial: ', '').replace('>', '').strip()
            name = "UNKNOWN"
            weight = info

            if polynomial_string == "TIMEOUT":
                name = "TIMEOUT"
            if polynomial_string == "Failed Projection":
                name = "Failed Projection"
            elif self.names_db is not None and polynomial_string in self.names_db['Polynomial'].values:
                name = self.names_db[self.names_db['Polynomial'] == polynomial_string]['Name'].iloc[0]

            histogram_data.append({"Polynomial": polynomial_string, "Frequency": weight, "Name": name})

        df = pd.DataFrame(histogram_data)
        df = df.sort_values(by='Frequency', ascending=False).reset_index(drop=True)
        # df['Polynomial'] = df['Polynomial'].apply(lambda x: x.replace('<Polynomial: ', '').replace('>', '').strip())
        if knotted_core:
            return df.iloc[0]['Name'], df.iloc[0]['Frequency'], df.iloc[0]['Polynomial']

        return df[['Name', 'Frequency', 'Polynomial']]
    

    def validate_arguments(self,matrix, max_nb_random_moves_III, max_nb_unsuccessfull_random_moves_III, timeout, flag_simplify_diagram):

        # Validate points
        if not isinstance(matrix, np.ndarray) or matrix.ndim != 2 or matrix.shape[1] != 3:
            raise ValueError("matrix must be a NumPy array of shape (n, 3).")
        
        # Validate max_moves_III and max_unsuccessful_moves_III
        if not isinstance(max_nb_random_moves_III, int) or max_nb_random_moves_III < 0:
            raise ValueError("max_moves_III must be a non-negative integer.")
        
        if not isinstance(max_nb_unsuccessfull_random_moves_III, int) or max_nb_unsuccessfull_random_moves_III < 0:
            raise ValueError("max_unsuccessful_moves_III must be a non-negative integer.")

        # Validate timeout
        if not isinstance(timeout, (int, float)) or timeout < 0:
            raise ValueError("timeout must be a non-negative number.")
        
        
        if not isinstance(flag_simplify_diagram, bool):
            raise ValueError("simplify_diagram must be a boolean value.")
        


    def analyze_curve(self, 
                              matrix, 
                              reduction_3d:bool=True, 
                              simplify_diagram:bool=True, 
                              max_nb_random_moves_III:int = 100_000,
                              max_nb_unsuccessfull_random_moves_III:int = 2_000,
                              timeout:int=1,
                              knotted_core:bool=False):

        self.validate_arguments(matrix, max_nb_random_moves_III, max_nb_unsuccessfull_random_moves_III, timeout, simplify_diagram)

        
        histogram= {}
        # Process each projection
        with self._suppress_cpp_output():
            for i, (dx, dy, dz) in enumerate(self.projections):
                weight = self.weights[i]

                polygontmp = kn.Polygon(matrix,self.flag_cyclic, self.flag_debug)
                polygontmp.set_closure(dx, dy, dz, self.closure_method)

                # Simplify 3D curve if necessary
                if reduction_3d:
                    if self.flag_debug:
                        print("Simplifying 3D curve...")
                    polygontmp.simplify_polygon(dx, dy, dz)
                    if self.flag_debug:
                        print(f"Polygon has {polygontmp.get_nb_points()} vertices.")


                # Project the polygon onto the plane
                flag_valid_projection = True

                try:
                    diagram = polygontmp.get_planar_diagram(dx, dy, dz, self.flag_planar)
                except Exception as e:
                    flag_valid_projection = False
                    print(f"Error: {str(e)}")

                if flag_valid_projection:
                    diagram.set_debug(self.flag_debug)
                    if simplify_diagram:
                        diagram.simplify()
                        if max_nb_random_moves_III > 0:
                            diagram.simplify_with_random_reidemeister_moves_III(max_nb_random_moves_III, max_nb_unsuccessfull_random_moves_III)
                            diagram.simplify()

                    # Evaluate polynomial invariants for the diagram
                    try:
                        jones = kn.PolynomialInvariant(diagram, self.flag_planar, self.flag_arrow_polynomial, self.flag_debug)
                        jones.set_timeout(timeout)
                        jones_polynomial = jones.get_polynomial_recursive()
                        self.update_histogram(histogram, str(jones_polynomial), weight)
                    except Exception as e:
                        print(f"Error: {str(e)}")
                        self.update_histogram(histogram, "TIMEOUT", weight)
                else:
                    self.update_histogram(histogram, "Failed Projection", weight)
        self.curve_analysis = self.create_histogram_dataframe(histogram,knotted_core)
        # return self.create_histogram_dataframe(histogram,knotted_core)
        return self.curve_analysis

    def get_polygon_slice(self,polygon_points, start_index, end_index, flag_cyclic):
        """
        Slices an nx3 numpy array representing a polygon based on the logic of the C++ get_polygon function.

        :param polygon_points: numpy array of shape (n, 3), each row representing a point (x, y, z).
        :param start_index: int, the start index for the slice.
        :param end_index: int, the end index for the slice.
        :param flag_cyclic: bool, indicates if the polygon is cyclic.
        :return: numpy array representing the sliced polygon.
        """
        n = polygon_points.shape[0]

        # Validating indices
        if flag_cyclic==False and start_index > end_index:
            raise ValueError("Start index must be less than or equal to end index for non-cyclic polygons.")
        if start_index < 0 or start_index >= n or end_index < 0 or (flag_cyclic==False and end_index >= n):
            raise ValueError("Invalid start index or end index.")

        # Adjusting indices for cyclic condition
        if end_index < start_index:
            end_index += n

        # Creating the slice
        indices = [i % n for i in range(start_index, end_index + 1)]
        return polygon_points[indices, :]

    def knotted_core(self,
                     matrix,
                     all_chains=False, 
                     closure:str=None, 
                     planar:bool=None, 
                     arrow_polynomial:bool=None, 
                     cyclic:bool=None, 
                     cyclic_input:bool=False,
                     kc_search_path:bool=False,
                     reduction_3d:bool=True,
                     simplify_diagram:bool=True,
                     timeout:int=1,
                     jones_method:str='recursive',
                     max_nb_random_moves_III:int = 100_000,
                     max_nb_unsuccessfull_random_moves_III:int = 2_000,
                     debug:bool=None
                     ):
        print('Performing knotted core analysis.')
        if all_chains:
            print('Sub-chain analysis enabled.')
        if kc_search_path:
            print('Search path output enabled.')

        if closure is None:
            closure = self.closure_method
        if planar is None:
            planar = self.flag_planar
        if arrow_polynomial is None:
            arrow_polynomial = self.flag_arrow_polynomial
        if cyclic is None:
            cyclic = self.flag_cyclic
        if debug is None:
            debug = self.flag_debug

        if cyclic and closure not in ('direct', 'rays'):
            raise ValueError(
                f"cyclic=True requires closure to be 'direct' or 'rays', got '{closure}'. "
                "A cyclic curve with 'open' closure aborts inside the C++ core."
            )

        json_data = self.names_db.to_json(orient='records')
        data_list = json.loads(json_data)
        formatted_str = ""
        for item in data_list:
            knot = item['Name']
            polynomial = item['Polynomial']
            formatted_str += f"{knot}\t{polynomial}\n"

        if type(self.projections) == np.ndarray:
            projections = self.projections.tolist()
        else:
            projections = self.projections
        if type(self.weights) == np.ndarray:
            weights = self.weights.tolist()
        else:
            weights = self.weights

        with self._suppress_cpp_output(debug):
            kc = kn.KnottedCore(matrix,
                                 formatted_str,
                                 projections,
                                 weights,
                                 all_chains,
                                 closure,
                                 planar,
                                 arrow_polynomial,
                                 cyclic,
                                 cyclic_input,
                                 kc_search_path,
                                 reduction_3d,
                                 simplify_diagram,
                                 timeout,
                                 jones_method,
                                 max_nb_random_moves_III,
                                 max_nb_unsuccessfull_random_moves_III,
                                 debug
                                 ).run()
        if all_chains:
            self.all_chains_analysis = pd.DataFrame.from_dict(kc[1])
            self.all_chains_analysis =self.all_chains_analysis.astype({'frequency':'float', 'index_first': 'int32', 'index_last': 'int32', 'length': 'int32'})
            if kc_search_path:
                # the C++ search-path map leaves the knot(oid)_type column empty,
                # so drop empty columns to keep the arrays the same length.
                search = {k: v for k, v in kc[2].items() if len(v) > 0}
                self.kc_search_path_analysis = pd.DataFrame(search)
        self.knotted_core_analysis = pd.DataFrame.from_dict(kc[0])
        return self.knotted_core_analysis, self.all_chains_analysis, self.kc_search_path_analysis

    # ------------------------------------------------------------------
    #  Plotting: translation of scripts/plot_knotted_core.R (ggplot2 -> plotnine)
    # ------------------------------------------------------------------
    # ColorBrewer "Set1" (qualitative), matching RColorBrewer used by the R script.
    _SET1 = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00',
             '#FFFF33', '#A65628', '#F781BF', '#999999']

    @staticmethod
    def _type_sort_key(t):
        """Natural sort of knot(oid) type names, matching plot_knotted_core.R."""
        m = re.match(r'([^0-9]*)([0-9]*)([^0-9]*)([0-9]*)(.*)', str(t))
        g1, g2, g3, g4, g5 = m.groups()
        return (g1, '*' in str(t), int(g2) if g2 else float('inf'),
                g3, int(g4) if g4 else float('inf'), g5)

    @staticmethod
    def _r_pretty(n=5):
        """Callable returning ~n 'nice' axis breaks, like R's scales::pretty_breaks."""
        def breaks(limits):
            lo, hi = min(limits), max(limits)
            rng = hi - lo
            if rng <= 0:
                return np.array([lo])
            mag = 10 ** np.floor(np.log10(rng / n))
            frac = (rng / n) / mag
            nice = 1 if frac < 1.5 else 2 if frac < 3 else 5 if frac < 7 else 10
            step = nice * mag
            start = np.floor(lo / step) * step
            stop = np.ceil(hi / step) * step
            return np.arange(start, stop + 0.5 * step, step)
        return breaks

    def _set1_colors(self, polynomials):
        """Assign Set1 colors to polynomials (in level order), ramping if >9,
        matching colorRampPalette(brewer.pal(9,'Set1')) in the R script."""
        from matplotlib.colors import LinearSegmentedColormap, to_hex
        n = len(polynomials)
        if n <= len(self._SET1):
            colors = self._SET1[:n]
        else:
            cmap = LinearSegmentedColormap.from_list('Set1', self._SET1)
            colors = [to_hex(cmap(i / (n - 1))) for i in range(n)]
        return dict(zip(polynomials, colors))

    def _prepare_plot_frame(self, df, N, no_transparency):
        """Rename C++/wrapper columns to the R data model and add 'middle'."""
        type_col = ('knot_type' if 'knot_type' in df.columns
                    else 'knotoid_type' if 'knotoid_type' in df.columns else None)
        rename = {'index_first': 'start_index', 'index_last': 'end_index'}
        if type_col:
            rename[type_col] = 'type'
        df = df.rename(columns=rename).copy()
        for c in ('start_index', 'end_index', 'length', 'frequency'):
            df[c] = pd.to_numeric(df[c])
        if 'type' not in df.columns:
            df['type'] = df['polynomial']
        if no_transparency:
            df['frequency'] = 1.0
        df['middle'] = (df['start_index'] + df['length'] / 2) % N
        return df

    def _fingerprint_boundaries(self, data, N):
        """White boundary segments between neighbouring cells of different type."""
        piv = (data.drop_duplicates(['start_index', 'end_index'])
                   .pivot(index='start_index', columns='end_index', values='polynomial'))
        starts, ends, M = piv.index.values, piv.columns.values, piv.values
        segs = []
        for i in range(len(starts) - 1):          # start vs start+1 -> vertical segment
            for j in range(len(ends)):
                a, b = M[i, j], M[i + 1, j]
                if isinstance(a, str) and isinstance(b, str) and a != b:
                    x = starts[i] + 0.5
                    segs.append((x, ends[j] - 0.5, x, ends[j] + 0.5))
        for i in range(len(starts)):              # end vs end+1 -> horizontal segment
            for j in range(len(ends) - 1):
                a, b = M[i, j], M[i, j + 1]
                if isinstance(a, str) and isinstance(b, str) and a != b:
                    y = ends[j] + 0.5
                    segs.append((starts[i] - 0.5, y, starts[i] + 0.5, y))
        return pd.DataFrame(segs, columns=['x', 'y', 'xend', 'yend'])

    def _search_path_line(self, path, N):
        """The blue search-path line (subchains of the initial type), split at wrap jumps."""
        first_poly = path['polynomial'].iloc[0]
        line = path[path['polynomial'] == first_poly].reset_index(drop=True)
        if len(line) < 2:
            return None
        ds = line['start_index'].diff().abs().fillna(0)
        de = line['end_index'].diff().abs().fillna(0)
        line = line.copy()
        line['grp'] = ((ds > N / 2) | (de > N / 2)).cumsum()
        return line

    def plot_knotted_core(self, output, cyclic=False, knotted_core=True,
                          search_path=False, boundaries=False, no_transparency=False):
        """Plot the subchain analysis, replicating scripts/plot_knotted_core.R.

        Requires knotted_core(all_chains=True) to have been run first (the
        search_path overlay additionally needs kc_search_path=True).

        :param output: output image filename.
        :param cyclic: disk matrix (polar: r=length, angle=middle) if True,
            otherwise fingerprint matrix (start index vs end index).
        :param knotted_core: overlay the knotted core(s) with yellow circles.
        :param search_path: overlay the search path (blue line, outlined cells).
        :param boundaries: draw white lines between domains of different type.
        :param no_transparency: ignore the dominant-type frequency (opaque tiles).
        """
        if self.all_chains_analysis is None:
            raise ValueError("No subchain analysis available. Run knotted_core(all_chains=True) first.")

        N = int(np.nanmax([
            pd.to_numeric(self.all_chains_analysis['index_first']).max() + 1,
            pd.to_numeric(self.all_chains_analysis['index_last']).max() + 1]))

        data = self._prepare_plot_frame(self.all_chains_analysis, N, no_transparency)

        # type factor ordered by natural sort => controls legend order and labels;
        # colors are assigned to polynomials in level order (as ggplot2 does).
        types = sorted(data['type'].astype(str).unique(), key=self._type_sort_key)
        data['type'] = pd.Categorical(data['type'].astype(str), categories=types, ordered=True)
        color_map = self._set1_colors(sorted(data['polynomial'].astype(str).unique()))
        legend = data[['type', 'polynomial']].drop_duplicates().sort_values('type')
        breaks = list(legend['polynomial'].astype(str))
        labels = [t if len(t) <= 40 else t[:40] + ' ...' for t in legend['type'].astype(str)]

        kc = None
        if knotted_core and self.knotted_core_analysis is not None:
            kc = self._prepare_plot_frame(self.knotted_core_analysis, N, no_transparency)

        path = None
        if search_path:
            if self.kc_search_path_analysis is None:
                raise ValueError("No search path available. Run knotted_core(all_chains=True, "
                                 "kc_search_path=True) first.")
            path = self._prepare_plot_frame(self.kc_search_path_analysis, N, no_transparency)

        if cyclic:
            fig = self._plot_disk(data, kc, path, boundaries, N, color_map, breaks, labels)
            fig.savefig(output, dpi=300, bbox_inches='tight')
            return fig

        p = self._plot_fingerprint(data, kc, path, boundaries, N, color_map, breaks, labels)
        p.save(output, width=10, height=6.667, dpi=300, verbose=False)
        return p

    def _plot_fingerprint(self, data, kc, path, boundaries, N, color_map, breaks, labels):
        import plotnine as p9
        p = (p9.ggplot(data, p9.aes(x='start_index', y='end_index'))
             + p9.geom_tile(p9.aes(fill='polynomial', alpha='frequency'))
             + p9.theme_bw()
             + p9.coord_fixed()
             + p9.xlab('start index') + p9.ylab('end index')
             + p9.scale_x_continuous(expand=(0, 0), limits=(-0.5, N + 0.5), breaks=self._r_pretty(5))
             + p9.scale_y_reverse(expand=(0, 0), limits=(N + 0.5, -0.5), breaks=self._r_pretty(5))
             + p9.scale_alpha(range=(0.1, 1), limits=(0, 1), guide=None)
             + p9.scale_fill_manual(values=color_map, breaks=breaks, labels=labels, name='')
             + p9.theme(legend_title=p9.element_blank()))
        if boundaries:
            seg = self._fingerprint_boundaries(data, N)
            if len(seg):
                p = p + p9.geom_segment(p9.aes(x='x', y='y', xend='xend', yend='yend'),
                                        data=seg, color='white', size=0.5, inherit_aes=False)
        if path is not None:
            p = p + p9.geom_tile(p9.aes(x='start_index', y='end_index'), data=path,
                                 color='black', alpha=0, size=0.2, inherit_aes=False)
            line = self._search_path_line(path, N)
            if line is not None:
                p = p + p9.geom_path(p9.aes(x='start_index', y='end_index', group='grp'),
                                     data=line, color='blue', size=0.5, inherit_aes=False)
        if kc is not None:
            p = p + p9.geom_point(p9.aes(x='start_index', y='end_index'), data=kc,
                                  color='yellow', fill='none', size=3, stroke=0.7, inherit_aes=False)
        return p

    def _disk_boundaries(self, data, N):
        """White boundary segments between neighbouring (start_index, length)
        cells of different type, in (middle, length) coordinates. Port of the
        cyclic branch of plot_knotted_core.R. Returns (m1, l1, m2, l2) tuples."""
        piv = (data.drop_duplicates(['start_index', 'length'])
                   .pivot(index='start_index', columns='length', values='polynomial'))
        starts, lengths, M = piv.index.to_numpy(), piv.columns.to_numpy(), piv.to_numpy()
        nr, nc = M.shape
        mid = lambda s, l: (s + l / 2.0) % N
        segs = []
        # start vs start+1 (cyclic), same length -> radial segment
        for i in range(nr):
            j = (i + 1) % nr
            for c in range(nc):
                a, b = M[i, c], M[j, c]
                if isinstance(a, str) and isinstance(b, str) and a != b:
                    m = mid(starts[i] + 0.5, lengths[c])
                    segs.append((m, lengths[c] - 0.5, m, lengths[c] + 0.5))
        # length vs length+1, same start -> arc segment
        for i in range(nr):
            for c in range(nc - 1):
                a, b = M[i, c], M[i, c + 1]
                if isinstance(a, str) and isinstance(b, str) and a != b:
                    m = mid(starts[i], lengths[c]); l = lengths[c] + 0.5
                    m1, m2 = m, m + 0.5
                    if m2 > N:
                        m1, m2 = m1 - N, m2 - N
                    segs.append((m1, l, m2, l))
        # length vs length+1 with start-1 -> arc segment
        for i in range(nr):
            ip = (i - 1) % nr
            for c in range(nc - 1):
                a, b = M[i, c], M[ip, c + 1]
                if isinstance(a, str) and isinstance(b, str) and a != b:
                    m = mid(starts[i], lengths[c]); l = lengths[c] + 0.5
                    m1, m2 = m - 0.5, m
                    if m1 < 0:
                        m1, m2 = m1 + N, m2 + N
                    segs.append((m1, l, m2, l))
        return segs

    def _plot_disk(self, data, kc, path, boundaries, N, color_map, breaks, labels):
        # plotnine has no coord_polar, so the disk (cyclic) matrix uses a
        # matplotlib polar axes styled to resemble the ggplot2 theme_bw output.
        from matplotlib.colors import to_rgba
        fig = plt.figure(figsize=(10, 6.667))
        ax = fig.add_subplot(111, projection='polar')
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)

        ang = (data['middle'].to_numpy() / N) * 2 * np.pi
        width = np.full(len(data), (1.0 / N) * 2 * np.pi)
        bottom = data['length'].to_numpy() - 0.5
        alphas = data['frequency'].clip(0.1, 1).to_numpy()
        facecolors = [to_rgba(color_map.get(str(p), '#999999'), a)
                      for p, a in zip(data['polynomial'], alphas)]
        ax.bar(ang, np.ones(len(data)), width=width, bottom=bottom,
               color=facecolors, edgecolor='none', align='center', linewidth=0)

        # boundaries between domains of different type (white)
        if boundaries:
            for m1, l1, m2, l2 in self._disk_boundaries(data, N):
                ax.plot([(m1 / N) * 2 * np.pi, (m2 / N) * 2 * np.pi], [l1, l2],
                        color='white', linewidth=0.6, solid_capstyle='round')

        # search path: outline the searched cells (black) and draw the path (blue)
        if path is not None:
            pang = (path['middle'].to_numpy() / N) * 2 * np.pi
            pwidth = np.full(len(path), (1.0 / N) * 2 * np.pi)
            pbottom = path['length'].to_numpy() - 0.5
            ax.bar(pang, np.ones(len(path)), width=pwidth, bottom=pbottom,
                   color='none', edgecolor='black', align='center', linewidth=0.3)
            line = self._search_path_line(path, N)
            if line is not None:
                for _, grp in line.groupby('grp'):
                    ax.plot((grp['middle'].to_numpy() / N) * 2 * np.pi,
                            grp['length'].to_numpy(), color='blue', linewidth=0.7)

        if kc is not None:
            kang = (kc['middle'].to_numpy() / N) * 2 * np.pi
            ax.plot(kang, kc['length'].to_numpy(), linestyle='none', marker='o',
                    markersize=8, markerfacecolor='none', markeredgecolor='yellow',
                    markeredgewidth=1.2)

        # angular ticks = middle point (0..N), theme_bw-like styling
        ticks = self._r_pretty(20)((0, N))
        ticks = ticks[(ticks >= 0) & (ticks < N)]
        ax.set_xticks((ticks / N) * 2 * np.pi)
        ax.set_xticklabels([f"{int(t)}" for t in ticks])
        ax.set_rlabel_position(112.5)   # keep radial (length) labels clear of the angular ones
        ax.set_facecolor('white')
        ax.grid(color='#b0b0b0', alpha=0.6, linewidth=0.5)
        ax.set_title('middle point (angle) / length (radius)', fontsize=10, pad=20)

        handles = [plt.Rectangle((0, 0), 1, 1, color=color_map[b]) for b in breaks]
        ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1.15, 0.5), frameon=False)
        fig.tight_layout()
        return fig

# %%
