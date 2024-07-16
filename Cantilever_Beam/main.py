# -*- coding: utf-8 -*-
"""
Author: Cássio Murakami
Project: Euler-Bernoulli Beam.
Title: main.py
"""
# Basic packages:
from config_packages import np, math, plt, cm, logging, data, module_name

# Configure the logging format and level (adjust as needed)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")

# Custom packages:
from src import mesh_generation, boundary_conditions, post_processing

def main():
    #0. Read data:
    logging.info(f"(0/6) Reading data from {module_name}.py ...")

    #1. Mesh generation:
    logging.info("(1/6) Generating mesh ...")
    mesh = mesh_generation.Mesh([data.x1, data.x2], data.h)

    #2. Assemble of global matrices and right hand side:
    logging.info("(2/6) Assembling global matrices and right hand side ...")

    node_dof = 2
    A = np.zeros([node_dof*mesh.ndof, node_dof*mesh.ndof])
    # Looping for each element:
    for element, local_element in enumerate(mesh.elements):
        element_nodes = mesh.elements_nodes[element]
        A_local = local_element.A_local()

        # Looping for each node of the element:
        for i, node_i in enumerate(element_nodes):
            for j, node_j in enumerate(element_nodes):
                global_i_base = node_dof * node_i
                local_i_base = node_dof * i
                global_j_base = node_dof * node_j
                local_j_base = node_dof * j
                
                # Looping for each degree of freedom of the node:
                for m, n in ((m, n) for m in range(node_dof) for n in range(node_dof)):
                    global_i = global_i_base + m
                    global_j = global_j_base + n
                    local_i = local_i_base + m
                    local_j = local_j_base + n
                
                    A[global_i][global_j] += A_local[local_i][local_j]

    F = np.zeros(node_dof*mesh.ndof)
    # Looping for each element:
    for element, local_element in enumerate(mesh.elements):
        element_nodes = mesh.elements_nodes[element]
        F_local = local_element.F_local()

        # Looping for each node of the element:
        for i, node_i in enumerate(element_nodes):
            global_i_base = node_dof * node_i
            local_i_base = node_dof * i

            # Looping for each degree of freedom of the node:
            for m in range(node_dof):
                global_i = global_i_base + m
                local_i = local_i_base + m

                F[global_i] += F_local[local_i] 


    # 3. Impose boundary conditions:
    logging.info("(3/6) Imposing boundary conditions ...")
    A, F, g = boundary_conditions.impose_boundary_conditions(A, F, mesh)

    # 4. Solve the algebraic problem:
    logging.info("(4/6) Solving the algebric problem ...")
    U = np.linalg.solve(A, F)
    # Lifting operation:
    U = U + g

    # Selecionar os itens de índice par usando list comprehension e enumerate
    y = np.array([item for index, item in enumerate(U) if index % 2 == 0])
    theta = -np.array([item for index, item in enumerate(U) if index % 2 == 1])

    # 5. Plotting the solution:
    logging.info("(5/6) Plotting the solution ...")
    # Plot the solution
    if data.plot_solution == 'y':
        fig = plt.figure(1)
        plt.plot(mesh.coord, data.u_analytical(mesh.coord), color = "black", label = "Analytical Solution")
        sc = plt.scatter(mesh.coord, y, c = y, cmap = 'jet', label = "Finite Element Method")
        plt.colorbar(sc)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend(loc="upper left")
        plt.title('Finite Element Solver: Euler-Bernoulli Beam')
        plt.draw()

        fig = plt.figure(2)
        plt.plot(mesh.coord, -data.grad_u_analytical(mesh.coord), color = "black", label = "Analytical Solution")
        sc = plt.scatter(mesh.coord, theta, c = theta, cmap = 'jet',label = "Finite Element Method")
        plt.colorbar(sc)
        plt.xlabel('x')
        plt.ylabel('Theta')
        plt.legend(loc="upper left")
        plt.title('Finite Element Solver: Euler-Bernoulli Beam')
        plt.draw()

    # 6. Computing the error:
    #logging.info("(6/6) Computing errors ...")
    #L2_error, H1_error = post_processing.compute_errors(U, mesh)
    #print(f"L2 Error = {L2_error}")
    #print(f"H1 Error = {H1_error}")

    # 7. Finish the program:
    plt.show()
    exit()

if __name__ == '__main__':
    main()