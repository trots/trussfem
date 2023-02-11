import math
import json
import vtk
import numpy as np
from typing import List, Dict, Optional


NULL_ID: int = -1


class Node:
    """
    A class used to represent a node of truss

    Attributes
    ----------
    id : int
        node ID
    x : float
        node X coordinate in metres
    y : float
        node Y coordinate in metres
    """

    id: int = NULL_ID
    x: float = 0.0
    y: float = 0.0


class Element:
    """
    A class used to represent an element of truss

    Attributes
    ----------
    id : int
        element ID
    profile_area: float
        element profile (cross-section) area in squared metres
    nodes: List[int]
        element's nodes IDs (constant size 2)
    """

    id: int = NULL_ID
    profile_area: float = 0.0
    nodes: List[int] = [NULL_ID, NULL_ID]


class Truss:
    """
    A class used to represent truss structure

    Attributes
    ----------
    nodes: List[Node]
        list of truss nodes (sorted by ID)
    elements: List[Element]
        list of truss elements (sorted by ID)
    """

    nodes: List[Node] = []
    elements: List[Element] = []


class Constraint:
    """
    A class used to represent a constraint in truss structure

    Attributes
    ----------
    id: int
        constraint ID
    node_id: int
        ID of a node where a constraint is applied
    contraint_x: bool
        constraint along X coordinate axis
    contraint_y: bool
        constraint along Y coordinate axis
    """

    id: int = NULL_ID
    node_id: int = NULL_ID
    constraint_x: bool = False
    constraint_y: bool = False


class Load:
    """
    A class used to represent a load in truss structure

    Attributes
    ----------
    id: int
        load ID
    node_id: int
        ID of a node where a load is applied
    load_x: float
        load in Newtons along X coordinate axis (positive if load direction coincides with X axis, negative otherwise)
    load_y: float
        load in Newtons along Y coordinate axis (positive if load direction coincides with Y axis, negative otherwise)
    """

    id: int = NULL_ID
    node_id: int = NULL_ID
    load_x: float = 0.0
    load_y: float = 0.0


class Result:
    """
    A class used to represent FEM results

    Attributes
    ----------
    displacements: np.ndarray
        displacements in metres of the each node of a truss
    nodal_forces: np.ndarray
        forces in Newtons of the each node of a truss
    strains: np.ndarray
        strains in metres of the each element of a truss
    stresses: np.ndarray
        stresses in Pascals of the each element of a truss
    """

    displacements: Optional[np.ndarray] = None
    nodal_forces: Optional[np.ndarray] = None
    strains: Optional[np.ndarray] = None
    stresses: Optional[np.ndarray] = None


class Task:
    """
    A class used to represent a task for this FEM implementation

    Attributes
    ----------
    name: str
        task name
    youngs_modulus: float
        Young's modulus of a truss material in Pascals
    poisson_coefficient: float
        poisson coefficient of a truss material (not used now)
    truss: Truss
        truss object
    truss: List[Constraint]
        list of constraints of a task (sorted by ID)
    loads: List[Load]
        list of loads of a task (sorted by ID)
    test_result: Result
        precise result of a task for test purposes (set as None if test not needed)
    """

    name: str = "Unknown task"
    youngs_modulus: float = 0.0
    poisson_coefficient: float = 0.0
    truss: Optional[Truss] = None
    constraints: List[Constraint] = []
    loads: List[Load] = []
    test_result: Optional[Result] = None


def decode_task(dct: Dict) -> Task:
    """Decodes JSON task to Task class object

    Parameters
    ----------
    dct : Dict
        JSON task description

    Returns
    -------
    task : Task
        Decoded task
    """

    task = Task()
    task.name = dct["name"]
    task.youngs_modulus = dct["youngs_modulus"]
    task.poisson_coefficient = dct["poisson_coefficient"]
    task.truss = Truss()

    for node in dct["truss"]["nodes"]:
        new_node = Node()
        new_node.id = node["id"]
        new_node.x = node["x"]
        new_node.y = node["y"]
        task.truss.nodes.append(new_node)

    for element in dct["truss"]["elements"]:
        new_element = Element()
        new_element.id = element["id"]
        new_element.nodes = [element["nodes"][0], element["nodes"][1]]
        new_element.profile_area = element["profile_area"]
        task.truss.elements.append(new_element)

    for constraint in dct["constraints"]:
        new_constraint = Constraint()
        new_constraint.id = constraint["id"]
        new_constraint.node_id = constraint["node_id"]
        new_constraint.constraint_x = constraint["constraint_x"]
        new_constraint.constraint_y = constraint["constraint_y"]
        task.constraints.append(new_constraint)

    for load in dct["loads"]:
        new_load = Load()
        new_load.id = load["id"]
        new_load.node_id = load["node_id"]
        new_load.load_x = load["load_x"]
        new_load.load_y = load["load_y"]
        task.loads.append(new_load)

    if "test_result" in dct:
        test = dct["test_result"]
        task.test_result = Result()
        task.test_result.displacements = test["displacements"]
        task.test_result.strains = test["strains"]
        task.test_result.stresses = test["stresses"]

    return task


def read_task(task_file_path: str) -> Task:
    """Reads a FEM task from file

    Parameters
    ----------
    task_file_path : str
        Path to a task file

    Returns
    -------
    task : Task
        Task object
    """

    with open(task_file_path, "r") as task_file:
        task_dict = json.load(task_file)
        task = decode_task(task_dict)
        return task


def element_length(element: Element, truss: Truss) -> float:
    """Returns element's length [6, p14]."""

    x_len = truss.nodes[element.nodes[1]].x - truss.nodes[element.nodes[0]].x
    y_len = truss.nodes[element.nodes[1]].y - truss.nodes[element.nodes[0]].y
    length = math.sqrt(x_len * x_len + y_len * y_len)
    return length


def build_local_d_matrices(truss: Truss) -> List[np.ndarray]:
    """ Creates and returns a list of D matrices for an each element of a truss. Naming as in [7, p11].
     Consider as B matrix in [5, p128]."""

    local_d_matrices = []

    for element in truss.elements:
        # Create B matrix
        length = element_length(element, truss)
        b_matrix = np.array([[1.0, 0.0], [1.0, length]])  # [5, p127]
        b_matrix_inv = np.linalg.inv(b_matrix)

        # Extend B matrix to (4,4) size for Y coordinates support
        b_matrix_inv = np.insert(b_matrix_inv, 1, [0.0], axis=0)
        b_matrix_inv = np.insert(b_matrix_inv, 3, [0.0], axis=0)
        b_matrix_inv = np.insert(b_matrix_inv, 1, [0.0], axis=1)
        b_matrix_inv = np.insert(b_matrix_inv, 3, [0.0], axis=1)

        # Create D matrix [5, p128]
        c_matrix = np.array([0, 0, 1, 0])
        local_d_matrix = c_matrix @ b_matrix_inv
        local_d_matrices.append(local_d_matrix)

    return local_d_matrices


def build_local_k_matrices(d_matrices: List[np.ndarray], youngs_modulus: float, truss: Truss) -> List[np.ndarray]:
    """ Creates a list of stiffness matrices K for an each element in a truss. K matrices are in local element's
    coordinates. K matrix size corresponds to [2, p12]. See also [5, p129], [7, p15]."""

    local_k_matrices: List[np.ndarray] = []

    for d_matrix, element in zip(d_matrices, truss.elements):
        coefficient = youngs_modulus * element.profile_area
        local_k_matrix = np.zeros((4, 4))
        local_k_matrix[0][0] = coefficient * d_matrix[2]
        local_k_matrix[2][2] = local_k_matrix[0][0]
        local_k_matrix[0][2] = coefficient * d_matrix[0]
        local_k_matrix[2][0] = local_k_matrix[0][2]
        local_k_matrices.append(local_k_matrix)

    return local_k_matrices


def build_c_matrices(truss: Truss) -> List[np.ndarray]:
    """ Creates a list of cosine matrices C for an each element of a truss. C matrix is used to convert K matrix to
     global coordinates. See [2, p22]."""

    c_matrices = []

    for element in truss.elements:
        global_x_vec = [1.0, 0.0]
        local_x_vec = [truss.nodes[element.nodes[1]].x - truss.nodes[element.nodes[0]].x,
                       truss.nodes[element.nodes[1]].y - truss.nodes[element.nodes[0]].y]
        cos = np.dot(local_x_vec, global_x_vec) / (np.linalg.norm(local_x_vec) * np.linalg.norm(global_x_vec))
        sin = math.sqrt(1 - cos * cos)
        c_matrix: np.ndarray = np.array([[cos, sin, 0.0, 0.0],
                                         [-sin, cos, 0.0, 0.0],
                                         [0.0, 0.0, cos, sin],
                                         [0.0, 0.0, -sin, cos]])
        c_matrices.append(c_matrix)

    return c_matrices


def convert_to_global_k_matrices(local_k_matrices: List[np.ndarray], c_matrices: List[np.ndarray]) -> List[np.ndarray]:
    """ Creates a list of stiffness matrices K in global coordinates system [2, p22]. """

    global_k_matrices: List[np.ndarray] = []

    for c_matrix, local_k_matrix in zip(c_matrices, local_k_matrices):
        global_k_matrix = c_matrix.T @ local_k_matrix @ c_matrix
        global_k_matrices.append(global_k_matrix)

    return global_k_matrices


def build_global_k_matrix(global_k_matrices: List[np.ndarray]) -> np.ndarray:
    """ Creates a one stiffness matrix in global coordinates system for the all elements [2, p26, p29]"""

    global_k_matrix_size = 4 * len(global_k_matrices)
    global_k_matrix = np.zeros((global_k_matrix_size, global_k_matrix_size))
    index = 0

    for k_matrix in global_k_matrices:
        lower = index * 4
        upper = index * 4 + 4
        global_k_matrix[lower:upper, lower:upper] = k_matrix
        index = index + 1

    return global_k_matrix


def bind_global_k_matrix(global_k_matrix: np.ndarray, truss: Truss) -> np.ndarray:
    """ Creates bindings between elements in a global stiffness matrix [2, p26, p30-31]. """

    h_matrix_rows = len(truss.elements) * 4
    h_matrix_cols = len(truss.nodes) * 2
    h_matrix: np.ndarray = np.zeros((h_matrix_rows, h_matrix_cols))
    e_matrix: np.ndarray = np.identity(2)
    element_index = 0

    for element in truss.elements:
        start_row = element_index * 4
        start_columns = [element.nodes[0] * 2, element.nodes[1] * 2]
        h_matrix[start_row:(start_row + 2), start_columns[0]:(start_columns[0] + 2)] = e_matrix
        h_matrix[(start_row + 2):(start_row + 4), start_columns[1]:(start_columns[1] + 2)] = e_matrix
        element_index = element_index + 1

    global_k_matrix = h_matrix.T @ global_k_matrix @ h_matrix
    return global_k_matrix


def apply_constraints(global_k_matrix: np.ndarray, constraints: List[Constraint]):
    """ Set constraints to a global K matrix [8, p19], [2, p32]"""

    for constraint in constraints:
        index_to_constraint = constraint.node_id * 2

        if constraint.constraint_x:
            global_k_matrix[index_to_constraint][index_to_constraint] = 1

            for i in range(global_k_matrix.shape[0]):
                if i == index_to_constraint:
                    continue
                global_k_matrix[index_to_constraint][i] = 0
                global_k_matrix[i][index_to_constraint] = 0

        index_to_constraint = index_to_constraint + 1

        if constraint.constraint_y:
            global_k_matrix[index_to_constraint][index_to_constraint] = 1

            for i in range(global_k_matrix.shape[0]):
                if i == index_to_constraint:
                    continue
                global_k_matrix[index_to_constraint][i] = 0
                global_k_matrix[i][index_to_constraint] = 0


def build_loads_vector(loads: List[Load], truss_nodes_count: int) -> np.ndarray:
    """ Creates loads vector with X and Y components for an each node. """
    loads_vector = np.zeros((truss_nodes_count * 2, 1))

    for load in loads:
        index = load.node_id * 2
        loads_vector[index] = load.load_x
        loads_vector[index + 1] = load.load_y

    return loads_vector


def solve(global_k_matrix: np.ndarray, loads_vector: np.ndarray) -> np.ndarray:
    """ Solve linear matrix equation and returns a nodal displacements vector. """

    displacements = np.linalg.lstsq(global_k_matrix, loads_vector, rcond=None)
    return displacements[0]


def get_global_element_displacements(element: Element, displacements: np.ndarray) -> np.ndarray:
    """ Returns displacements vector for an element. """

    element_displacements = np.zeros((4, 1))
    node_0_row = element.nodes[0] * 2
    element_displacements[0:2] = displacements[[node_0_row, node_0_row + 1]]
    node_1_row = element.nodes[1] * 2
    element_displacements[2:4] = displacements[[node_1_row, node_1_row + 1]]
    return element_displacements


def get_local_element_displacements(element: Element, c_matrix: np.ndarray, displacements: np.ndarray) -> np.ndarray:
    """ Convert a global element's displacements vector to local coordinates system [2, p34]. """

    element_displacements = get_global_element_displacements(element, displacements)
    element_displacements = c_matrix @ element_displacements
    return element_displacements


def build_nodal_forces(local_k_matrices: List[np.ndarray], c_matrices: List[np.ndarray], displacements: np.ndarray,
                       truss: Truss) -> List[np.ndarray]:
    """ Creates a nodal forces vector for an each element [2, p34]. """

    nodal_forces_list = []

    for k_matrix, c_matrix, element in zip(local_k_matrices, c_matrices, truss.elements):
        element_displacements = get_local_element_displacements(element, c_matrix, displacements)
        nodal_forces = k_matrix @ element_displacements
        nodal_forces_list.append(nodal_forces)

    return nodal_forces_list


def build_strains_vector(truss: Truss, d_matrices: List[np.ndarray], c_matrices: List[np.ndarray],
                         displacements: np.ndarray) -> np.ndarray:
    """ Creates a strains vector [7, p10]. """

    strains = []

    for element, d_matrix, c_matrix in zip(truss.elements, d_matrices, c_matrices):
        element_displacements = get_local_element_displacements(element, c_matrix, displacements)
        element_strain = d_matrix @ element_displacements
        strains.append(element_strain[0])

    return np.asarray(strains)


def build_stresses_vector(strains: np.ndarray, youngs_modulus: float) -> np.ndarray:
    """ Creates a stresses vector [7, p12]. """

    stresses = youngs_modulus * strains
    return stresses


def write_displacements_to_vtk(truss: Truss, displacements: np.ndarray, filename: str):
    points = vtk.vtkPoints()

    for node in truss.nodes:
        points.InsertNextPoint([node.x, node.y, 0.0])

    vtk_scalars_array_x = vtk.vtkDoubleArray()
    vtk_scalars_array_x.SetName("displacements_x")
    vtk_scalars_array_y = vtk.vtkDoubleArray()
    vtk_scalars_array_y.SetName("displacements_y")

    for i in range(0, displacements.shape[0], 2):
        vtk_scalars_array_x.InsertNextValue(displacements[i])
        vtk_scalars_array_y.InsertNextValue(displacements[i + 1])

    lines = vtk.vtkCellArray()

    for element in truss.elements:
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, element.nodes[0])
        line.GetPointIds().SetId(1, element.nodes[1])
        lines.InsertNextCell(line)

    lines_poly_data = vtk.vtkPolyData()
    lines_poly_data.SetPoints(points)
    lines_poly_data.SetLines(lines)
    lines_poly_data.GetPointData().AddArray(vtk_scalars_array_x)
    lines_poly_data.GetPointData().AddArray(vtk_scalars_array_y)

    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(lines_poly_data)
    writer.SetFileName(filename)
    writer.Update()


def write_elements_data_to_vtk(truss: Truss, data_array: np.ndarray, data_name: str, data_filename: str):
    points = vtk.vtkPoints()

    for node in truss.nodes:
        points.InsertNextPoint([node.x, node.y, 0.0])

    vtk_data_array = vtk.vtkDoubleArray()
    vtk_data_array.SetName(data_name)

    for data_item in data_array:
        vtk_data_array.InsertNextValue(data_item)

    lines = vtk.vtkCellArray()

    for element in truss.elements:
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, element.nodes[0])
        line.GetPointIds().SetId(1, element.nodes[1])
        lines.InsertNextCell(line)

    lines_poly_data = vtk.vtkPolyData()
    lines_poly_data.SetPoints(points)
    lines_poly_data.SetLines(lines)
    lines_poly_data.GetCellData().AddArray(vtk_data_array)

    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(lines_poly_data)
    writer.SetFileName(data_filename)
    writer.Update()
