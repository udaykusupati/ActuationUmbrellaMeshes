import umbrella_mesh
from bending_validation import suppress_stdout as so
from load_jsondata import read_data
import numpy as np
import py_newton_optimizer
from equilibrium_solve_analysis import EquilibriumSolveAnalysis
from pipeline_helper import UmbrellaOptimizationCallback, allEnergies, allGradientNorms, allDesignObjectives, allDesignGradientNorms, set_joint_vector_field, show_center_joint_normal, show_joint_normal

    
def parse_input(input_path, handleBoundary = False, resolution = 10, handlePivots = True, isHex = False, material_params = None, circular_cs = False, use_target_surface = True):
    if material_params is not None:
        input_data, io = read_data(filepath = input_path, handleBoundary = handleBoundary, handlePivots = handlePivots, isHex = isHex, material_params=material_params, circularCS = circular_cs, useTargetSurface = use_target_surface)
    else:
        input_data, io = read_data(filepath = input_path, handleBoundary = handleBoundary, handlePivots = handlePivots, isHex = isHex, circularCS = circular_cs, useTargetSurface = use_target_surface)
    import mesh
    target_mesh = None
    if len(input_data['target_v']) != 0:
        target_mesh = mesh.Mesh(input_data['target_v'], input_data['target_f'])
    
    # for vid in range(len(input_data['vertices'])):
        # if input_data['v_labels'][vid] == "TH":
        #     print("TH", input_data['vertices'][vid])
        # if input_data['v_labels'][vid] == "BH":
        #     print("BH", input_data['vertices'][vid])
    # print(input_data['vertices'][3], input_data['vertices'][20])
    # print(input_data['v_labels'][3], input_data['v_labels'][20])
    # print(input_data['e_labels'][30])
    # for seg_id in range(len(io.segments)):
    #     print(seg_id,io.segments[seg_id].type, io.segments[seg_id].normal )
    curr_um = umbrella_mesh.UmbrellaMesh(io, resolution)
    plate_thickness = io.material_params[-2]
    target_height_multiplier = input_data['target_spacing_factor']
    return io, input_data, target_mesh, curr_um, plate_thickness, target_height_multiplier

def configure_umbrella_pre_deployment(curr_um, thickness, target_height_multiplier):
    curr_um.deploymentForceType = umbrella_mesh.DeploymentForceType.LinearActuator
    if isinstance(target_height_multiplier, np.ndarray):
        curr_um.targetDeploymentHeightVector = thickness * target_height_multiplier
    else:
        curr_um.targetDeploymentHeight = thickness * target_height_multiplier
    curr_um.repulsionEnergyWeight = 0
    curr_um.attractionWeight = 1e-1
    curr_um.setHoldClosestPointsFixed(False)
    curr_um.scaleInputPosWeights(0.99999)

    curr_um.angleBoundEnforcement = umbrella_mesh.AngleBoundEnforcement.Penalty
    # curr_um.angleBoundEnforcement = umbrella_mesh.AngleBoundEnforcement.Hard

def break_input_angle_symmetry(curr_um):
    dof = curr_um.getDoFs()
    for i in range(curr_um.numJoints()):
        # if (curr_um.joint(i).jointType() == umbrella_mesh.JointType.X):
        dof[curr_um.dofOffsetForJoint(i) + 6] = 1e-3
    curr_um.setDoFs(dof)

def insert_randomness(curr_um, zPerturbationEpsilon = 1e-4):
    dof = np.array(curr_um.getDoFs())
    zCoordDoFs = np.array(curr_um.jointPositionDoFIndices())[2::3]
    dof[zCoordDoFs] += 2 * zPerturbationEpsilon * (np.random.random_sample(len(zCoordDoFs)) - 0.5)
    curr_um.setDoFs(dof)

'''
def staged_deployment(curr_um, weights, eqm_callback, OPTS, fixedVars, elasticEnergyIncreaseFactorLimit = 2.5):
    for weight in weights:
        curr_um.uniformDeploymentEnergyWeight = weight
        with so(): 
            results = umbrella_mesh.compute_equilibrium(curr_um, callback = eqm_callback, options = OPTS, fixedVars = fixedVars, elasticEnergyIncreaseFactorLimit=elasticEnergyIncreaseFactorLimit)
    
    curr_um.angleBoundEnforcement = umbrella_mesh.AngleBoundEnforcement.Hard
    with so(): 
        results = umbrella_mesh.compute_equilibrium(curr_um, callback = eqm_callback, options = OPTS, fixedVars = fixedVars, elasticEnergyIncreaseFactorLimit=elasticEnergyIncreaseFactorLimit)
    return results
'''

def configure_umbrella_optimization(curr_um, bdryMultiplier = 1.0):
    # ### Initialize Design Optimization
    curr_um.uniformDeploymentEnergyWeight = 1e0
    curr_um.repulsionEnergyWeight = 0
    curr_um.attractionWeight = 1e-3
    curr_um.setHoldClosestPointsFixed(True) # Don't let closest point's drift away from reasonable values with the low attraction weight.
    curr_um.scaleInputPosWeights(0.1, bdryMultiplier)


def configure_umbrella_true_equlibrium(curr_um, thickness, target_height_multiplier):
    curr_um.uniformDeploymentEnergyWeight = 1e0
    curr_um.targetDeploymentHeight = thickness * target_height_multiplier 
    curr_um.attractionWeight = 0

def configure_umbrella_undeployment_step_one(curr_um, thickness, target_height_multiplier, undeployment_multiplier = 10):
    curr_um.uniformDeploymentEnergyWeight = 1e0
    curr_um.targetDeploymentHeight = thickness * target_height_multiplier * undeployment_multiplier
    curr_um.attractionWeight = 0

def configure_umbrella_undeployment_step_two(curr_um):
    curr_um.uniformDeploymentEnergyWeight = 0
    curr_um.attractionWeight = 0

def configure_design_optimization_umbrella(uo):
    uo.equilibriumOptimizer.options.verbose = 1
    #uo.equilibriumOptimizer.options.verboseWorkingSet = True
    uo.equilibriumOptimizer.options.gradTol = 1e-10
    # Hold the closest points fixed in the target-attraction term of the equilibrium solve:
    # this seems to make the design optimization much more robust.
    uo.setHoldClosestPointsFixed(True, False)

'''    
def deploy_umbrella_pin_rigid_motion(curr_um, plate_thickness, target_height_multiplier, view = None, colors = None, releaseActuation = False, analysis = False, dep_weights = np.logspace(-3, 0, 4)):
    use_pin = True
    driver = curr_um.centralJoint()
    driverj = curr_um.joint(driver)
    jdo = curr_um.dofOffsetForJoint(driver)
    fixedVars = (list(range(jdo, jdo + 6)) if use_pin else []) + curr_um.rigidJointAngleDoFIndices()
    
    OPTS = py_newton_optimizer.NewtonOptimizerOptions()
    OPTS.gradTol = 1e-8
    OPTS.verbose = 1
    OPTS.beta = 1e-6
    OPTS.niter = 300
    OPTS.verboseNonPosDef = False
    
    eqays = EquilibriumSolveAnalysis(curr_um)
    def eqm_callback(prob, i):
        if analysis: eqays.record(prob)
        if (i % 1 == 0):
            if view is not None:
                view.update(scalarField = colors)
                view.show()


    configure_umbrella_pre_deployment(curr_um, plate_thickness, target_height_multiplier)
    if curr_um.getTargetSurface() is None:
        curr_um.attractionWeight = 0
    
    break_input_angle_symmetry(curr_um)
    if releaseActuation:
        results = staged_deployment(curr_um, [0], eqm_callback, OPTS, fixedVars) ## Use more/slower steps if the deployment is hard/tangled
    else:
        results = staged_deployment(curr_um, dep_weights, eqm_callback, OPTS, fixedVars) ## Use more/slower steps if the deployment is hard/tangled
        # results = staged_deployment(curr_um, np.logspace(-5, 0, 6), eqm_callback, OPTS, fixedVars) ## Use more/slower steps if the deployment is hard/tangled
    return results.success, eqays
'''
    
def get_material_params(mat):
    E, nu = None, None
    if mat == "PP":
        E, nu = 1400, 0.35
    elif mat == "POM-C":
        E, nu = 3000, 0.44
    else:
        assert 0, "Unknown material"
    return E, nu



# ==============================

def staged_deployment(curr_um, weights, eqm_callback, OPTS, fixedVars, elasticEnergyIncreaseFactorLimit = 2.5):
    for weight in weights:
        if isinstance(weight, np.ndarray):
            curr_um.deploymentEnergyWeight = weight
        else:
            curr_um.uniformDeploymentEnergyWeight = weight

        with so(): 
            results = umbrella_mesh.compute_equilibrium(curr_um, callback = eqm_callback, options = OPTS, fixedVars = fixedVars, elasticEnergyIncreaseFactorLimit=elasticEnergyIncreaseFactorLimit)
    
    curr_um.angleBoundEnforcement = umbrella_mesh.AngleBoundEnforcement.Hard
    with so(): 
        results = umbrella_mesh.compute_equilibrium(curr_um, callback = eqm_callback, options = OPTS, fixedVars = fixedVars, elasticEnergyIncreaseFactorLimit=elasticEnergyIncreaseFactorLimit)
    return results


def deploy_umbrella_pin_rigid_motion(curr_um, plate_thickness, target_height_multiplier, view = None, colors = None, releaseActuation = False, analysis = False, dep_weights = np.logspace(-3, 0, 4)):
    use_pin = True
    driver = curr_um.centralJoint()
    jdo = curr_um.dofOffsetForJoint(driver)
    fixedVars = (list(range(jdo, jdo + 6)) if use_pin else []) + curr_um.rigidJointAngleDoFIndices()
    
    OPTS = py_newton_optimizer.NewtonOptimizerOptions()
    OPTS.gradTol = 1e-8
    OPTS.verbose = 1
    OPTS.beta = 1e-6
    OPTS.niter = 600 # default value is 300
    OPTS.verboseNonPosDef = False
    
    eqays = EquilibriumSolveAnalysis(curr_um)
    def eqm_callback(prob, i):
        if analysis:
            eqays.record(prob)
        if view is not None:
            view.update(scalarField = colors)
            view.show()


    configure_umbrella_pre_deployment(curr_um, plate_thickness, target_height_multiplier)
    if curr_um.getTargetSurface() is None:
        curr_um.attractionWeight = 0
        # [RK] play with attractionWeight
    
    break_input_angle_symmetry(curr_um)
    if releaseActuation:
        results = staged_deployment(curr_um, [0], eqm_callback, OPTS, fixedVars) ## Use more/slower steps if the deployment is hard/tangled
    else:
        results = staged_deployment(curr_um, dep_weights, eqm_callback, OPTS, fixedVars) ## Use more/slower steps if the deployment is hard/tangled
        # results = staged_deployment(curr_um, np.logspace(-5, 0, 6), eqm_callback, OPTS, fixedVars) ## Use more/slower steps if the deployment is hard/tangled
    return results.success, eqays
