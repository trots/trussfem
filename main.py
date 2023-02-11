import sys
import numpy as np
import trussfem as fem

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Not enough arguments. Specify the task file")
        exit(0)

    result: fem.Result = fem.Result()

    task_file_path = sys.argv[1]
    print("Read input task:", task_file_path)
    task = fem.read_task(task_file_path)
    print("Task:", task.name)
    print("Nodes:", len(task.truss.nodes))
    print("Elements:", len(task.truss.elements))

    print("Build local D matrices...", end=" ")
    d_matrices = fem.build_local_d_matrices(task.truss)
    print("Ok")

    print("Build local K matrices...", end=" ")
    local_k_matrices = fem.build_local_k_matrices(d_matrices, task.youngs_modulus, task.truss)
    print("Ok")

    print("Build C matrices...", end=" ")
    c_matrices = fem.build_c_matrices(task.truss)
    print("Ok")

    print("Convert to global K matrices...", end=" ")
    global_k_matrices = fem.convert_to_global_k_matrices(local_k_matrices, c_matrices)
    print("Ok")

    print("Build global K matrix...", end=" ")
    global_k_matrix = fem.build_global_k_matrix(global_k_matrices)
    print("Ok")

    print("Bind global K matrix...", end=" ")
    global_k_matrix = fem.bind_global_k_matrix(global_k_matrix, task.truss)
    print("Ok")

    print("Apply constraints...", end=" ")
    fem.apply_constraints(global_k_matrix, task.constraints)
    print("Ok")

    print("Build loads vector...", end=" ")
    loads_vector = fem.build_loads_vector(task.loads, len(task.truss.nodes))
    print("Ok")

    print("Solve...", end=" ")
    result.displacements = fem.solve(global_k_matrix, loads_vector)
    print("Ok")

    print("Write displacements to vtk file...", end=" ")
    fem.write_displacements_to_vtk(task.truss, result.displacements, "displacements.vtk")
    print("Ok")

    print("Build nodal forces...", end=" ")
    result.nodal_forces = fem.build_nodal_forces(local_k_matrices, c_matrices, result.displacements, task.truss)
    print("Ok")

    print("Build strains...", end=" ")
    result.strains = fem.build_strains_vector(task.truss, d_matrices, c_matrices, result.displacements)
    print("Ok")

    print("Write strains to vtk file...", end=" ")
    fem.write_elements_data_to_vtk(task.truss, result.strains, "strains", "strains.vtk")
    print("Ok")

    print("Build stresses...", end=" ")
    result.stresses = fem.build_stresses_vector(result.strains, task.youngs_modulus)
    print("Ok")

    print("Write stresses to vtk file...", end=" ")
    fem.write_elements_data_to_vtk(task.truss, result.stresses, "stresses", "stresses.vtk")
    print("Ok")

    if task.test_result is not None:
        print("Result testing...")
        all_tests_ok = True
        test_ok = np.array_equiv(result.displacements, task.test_result.displacements)

        if not test_ok:
            print("ERROR: Displacements test failed")
            all_tests_ok = False

        test_ok = np.array_equiv(result.strains, task.test_result.strains)

        if not test_ok:
            print("ERROR: Strains test failed")
            all_tests_ok = False

        test_ok = np.array_equiv(result.stresses, task.test_result.stresses)

        if not test_ok:
            print("ERROR: Stresses test failed")
            all_tests_ok = False

        if all_tests_ok:
            print("Test succeed")
        else:
            print("Test failed")

    print("Completed")
