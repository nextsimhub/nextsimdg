"""GDB pretty-plotting commands for neXtSIM-DG field types.

Provides following GDB commands:

- ``display_marray``:   plot a component of a ``ModelArray`` (H, DG, VERTEX, ...).
- ``display_dgvector``: plot a component of a ``DGVector<N>`` (Eigen matrix).
- ``display_cgvector``: plot a ``CGVector<degree>`` (flat CG node array).

``display_marray``/``display_dgvector`` take ``-c/--component`` (default 0);
``display_dgvector``/``display_cgvector`` require ``--nx``/``--ny``, and
``display_cgvector`` additionally requires ``--degree``.
"""

import argparse

import matplotlib.pyplot as plt
import mplcursors
import numpy as np

import gdb


def cell_field_plot(data, component=0):
    """Plot one component of a cell-centred field (e.g. H, DG or DGVector)."""
    fig, ax = plt.subplots()

    field = data[:, :, component]
    ax.pcolor(field, edgecolors="k", shading="flat")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_title(f"Component {component}")

    # mplcursors can't pick a pcolor mesh, so use invisible points at cell centres
    x, y = np.meshgrid(np.arange(field.shape[1]) + 0.5, np.arange(field.shape[0]) + 0.5)
    sc = ax.scatter(x, y, c=field, s=100)

    cursor = mplcursors.cursor(sc, hover=True)

    @cursor.connect("add")
    def _(sel):
        i, j = np.unravel_index(sel.index, field.shape)
        sel.annotation.set_text(f"({j}, {i}) = {field[i, j]:g}")

    fig.tight_layout()
    plt.show()
    plt.close(fig)


def vertex_field_plot(data, component=0):
    """Plot the vertex field as coloured points."""
    fig, ax = plt.subplots()

    field = data[:, :, component]

    x, y = np.meshgrid(np.arange(field.shape[1]), np.arange(field.shape[0]))

    # Single component, used for colouring
    sc = ax.scatter(x, y, c=field, cmap="viridis", s=100)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_title(f"Component {component}")

    cursor = mplcursors.cursor(sc, hover=True)

    @cursor.connect("add")
    def _(sel):
        i, j = np.unravel_index(sel.index, field.shape)
        sel.annotation.set_text(f"({j}, {i}) = {field[i, j]:g}")

    fig.tight_layout()
    plt.show()
    plt.close(fig)


def cg_field_plot(data, nx, ny, degree):
    """Plot a CG field, colouring every CG node and drawing the cell grid.

    data is the flat array of CG nodes of shape ((degree*ny+1) * (degree*nx+1),).
    The axes enumerate cell indices: integer coordinates are cell boundaries,
    so cell (i, j) spans x in [i, i+1], y in [j, j+1]; interior CG nodes fall
    between the grid lines.
    """
    fig, ax = plt.subplots()

    n_nodes_x = degree * nx + 1
    n_nodes_y = degree * ny + 1
    field = data.reshape((n_nodes_y, n_nodes_x))  # row = y (j), col = x (i)

    # Node position in CELL units: x = i/degree, y = j/degree (x is the fast index)
    x, y = np.meshgrid(np.arange(n_nodes_x) / degree, np.arange(n_nodes_y) / degree)

    # Colour every CG node by its value (includes edge + interior nodes)
    sc = ax.scatter(x, y, c=field, cmap="viridis", s=60)

    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)
    ax.set_xticks(np.arange(nx + 1))  # integer ticks -> cell boundaries
    ax.set_yticks(np.arange(ny + 1))
    ax.grid(True, color="k", linewidth=0.7)  # grid lines = cell boundaries
    ax.set_xlabel("cell index x")
    ax.set_ylabel("cell index y")
    ax.set_aspect("equal")
    ax.set_title(f"CG field (degree {degree})")

    cursor = mplcursors.cursor(sc, hover=True)

    @cursor.connect("add")
    def _(sel):
        j, i = np.unravel_index(sel.index, field.shape)  # field is (y, x)
        # Node position in cell units
        cx, cy = i / degree, j / degree
        is_vertex = (i % degree == 0) and (j % degree == 0)
        kind = (
            "vertex"
            if is_vertex
            else ("edge" if (i % degree == 0) or (j % degree == 0) else "interior")
        )
        sel.annotation.set_text(
            f"node ({i},{j}) [{kind}]\ncell ({cx:.3g}, {cy:.3g}) = {field[j, i]:g}"
        )

    fig.tight_layout()
    plt.show()
    plt.close(fig)


def read_to_list(ptr, n0, n1):
    """Read the Eigen buffer into a Python list of lists.

    In nextsim we use Eigen in a row-major order.
    """
    # Read the memory buffer
    size = n0 * n1
    data = []
    for i in range(size):
        data.append(float(ptr[i]))
    return data


def dispatch_plot(array_type):
    """Return the plotting function appropriate for a ModelArray type."""
    if array_type in (
        gdb.parse_and_eval("ModelArray::Type::H"),
        gdb.parse_and_eval("ModelArray::Type::DG"),
    ):
        return cell_field_plot
    if array_type == gdb.parse_and_eval("ModelArray::Type::VERTEX"):
        return vertex_field_plot
    msg = f"Unsupported array type: {array_type}"
    raise ValueError(msg)


class DisplayMArray(gdb.Command):
    """Make a plot of the data in a ModelArray."""

    def __init__(self):
        super().__init__("display_marray", gdb.COMMAND_USER)

    @staticmethod
    def _make_argument_parser():
        """Build the argument parser for the display_marray command."""
        parser = argparse.ArgumentParser(prog="display_marray")
        parser.add_argument("variable", help="The ModelArray variable to plot")
        parser.add_argument(
            "-c",
            "--component",
            type=int,
            default=0,
            help="Index of the field component to plot (default: 0)",
        )
        return parser

    @staticmethod
    def _parse_arguments(arg):
        """Parse the arguments of the display_marray GDB command.

        Only gdb.parse_and_eval is allowed to throw errors since we can run anything
        through it, which may result in various exceptions being thrown. Therefore we
        try and parse as much as possible with argparse to make sure the inputs are
        sensible.
        """
        args = DisplayMArray._make_argument_parser().parse_args(gdb.string_to_argv(arg))
        return args.variable, args.component

    def invoke(self, arg, from_tty):
        """Invoke the command."""
        # Parse the command line arguments
        variable, component = self._parse_arguments(arg)

        # Get the variable into a gdb.Value
        val = gdb.parse_and_eval(variable)

        marray_type = val["type"]

        # Very Hacky way to get the dimensions of the array
        nx = gdb.parse_and_eval(f"{variable}.dimensions()[0]").referenced_value()
        ny = gdb.parse_and_eval(f"{variable}.dimensions()[1]").referenced_value()

        # Get the info about the underlying Eigen buffer
        storage_info = val["m_data"]["m_storage"]

        ptr = storage_info["m_data"]
        row_size = storage_info["m_cols"]
        col_size = storage_info["m_rows"]  # This is also the number of components

        # Establish a numpy array on the data
        data = np.array(read_to_list(ptr, col_size, row_size))

        # We make assumptions here on the order
        data = data.reshape((ny, nx, row_size))

        dispatch_plot(marray_type)(data, component)


class DisplayDGVector(gdb.Command):
    """Make a plot of the data in a DGVector<N>.

    A DGVector is a plain Eigen matrix (rows = cells, cols = DG components)
    without ModelArray metadata, so the grid dimensions are supplied by the
    user via --nx/--ny (GDB expressions) and the data layout is derived from
    the matrix's rows()/cols().
    """

    def __init__(self):
        super().__init__("display_dgvector", gdb.COMMAND_USER)

    @staticmethod
    def _make_argument_parser():
        """Build the argument parser for the display_dgvector command."""
        parser = argparse.ArgumentParser(prog="display_dgvector")
        parser.add_argument("variable", help="The DGVector variable to plot")
        parser.add_argument(
            "-c",
            "--component",
            type=int,
            default=0,
            help="Index of the DG component to plot (default: 0)",
        )
        parser.add_argument(
            "--nx",
            type=str,
            required=True,
            help="Grid extent in x. A GDB expression, evaluated with parse_and_eval "
            "(e.g. '20', 'smesh.nx', or 'metadata.getLocalExtentX() + 2').",
        )
        parser.add_argument(
            "--ny",
            type=str,
            required=True,
            help="Grid extent in y. A GDB expression, evaluated with parse_and_eval.",
        )
        return parser

    @staticmethod
    def _parse_arguments(arg):
        """Parse the arguments of the display_dgvector GDB command.

        Only gdb.parse_and_eval is allowed to throw errors since we can run anything
        through it, which may result in various exceptions being thrown. Therefore we
        try and parse as much as possible with argparse to make sure the inputs are
        sensible.

        Returns (variable, component, nx_expr, ny_expr) where nx/ny are left as
        strings to be evaluated in invoke (in the current frame).
        """
        args = DisplayDGVector._make_argument_parser().parse_args(
            gdb.string_to_argv(arg)
        )
        return args.variable, args.component, args.nx, args.ny

    def invoke(self, arg, from_tty):
        """Invoke the command."""
        # Parse the command line arguments
        variable, component, nx_expr, ny_expr = self._parse_arguments(arg)

        # Get the variable into a gdb.Value
        val = gdb.parse_and_eval(variable)

        # DGVector dimensions live behind rows()/cols() methods (not fields),
        # so evaluate them through GDB in the current frame.
        n_cells = int(gdb.parse_and_eval(f"{variable}.rows()"))
        n_comps = int(gdb.parse_and_eval(f"{variable}.cols()"))

        # nx/ny are GDB expressions supplied by the user
        nx = int(gdb.parse_and_eval(nx_expr))
        ny = int(gdb.parse_and_eval(ny_expr))

        if nx * ny != n_cells:
            msg = (
                f"nx*ny = {nx}*{ny} = {nx * ny} does not match the DGVector "
                f"cell count ({n_cells})"
            )
            raise ValueError(msg)

        # Sanity check the requested component
        if not 0 <= component < n_comps:
            msg = f"Component {component} is out of range for {n_comps} DG components"
            raise ValueError(msg)

        # Get the info about the underlying Eigen buffer (DGVector is the matrix)
        storage_info = val["m_storage"]
        ptr = storage_info["m_data"]

        # Establish a numpy array on the data
        data = np.array(read_to_list(ptr, n_cells, n_comps))

        # DGVector stores (cells x components) row-major; cells are (j, i)-ordered
        data = data.reshape((ny, nx, n_comps))

        cell_field_plot(data, component)


class DisplayCGVector(gdb.Command):
    """Make a plot of the data in a CGVector<degree>.

    A CGVector is a flat Eigen column vector of CG nodes of length
    (degree*nx+1)*(degree*ny+1), ordered with x fastest from the bottom-left.
    The grid extents and the CG degree are user-supplied (--nx/--ny as GDB
    expressions, --degree as a plain int); the layout is validated against
    size() read from GDB.
    """

    def __init__(self):
        super().__init__("display_cgvector", gdb.COMMAND_USER)

    @staticmethod
    def _make_argument_parser():
        """Build the argument parser for the display_cgvector command."""
        parser = argparse.ArgumentParser(prog="display_cgvector")
        parser.add_argument("variable", help="The CGVector variable to plot")
        parser.add_argument(
            "--nx",
            type=str,
            required=True,
            help="Grid extent in x (cell count). GDB expression via parse_and_eval.",
        )
        parser.add_argument(
            "--ny",
            type=str,
            required=True,
            help="Grid extent in y (cell count). GDB expression via parse_and_eval.",
        )
        parser.add_argument(
            "--degree",
            type=int,
            required=True,
            help="CG polynomial degree (number of sub-intervals per cell).",
        )
        return parser

    @staticmethod
    def _parse_arguments(arg):
        """Parse the arguments of the display_cgvector GDB command.

        Only gdb.parse_and_eval is allowed to throw errors since we can run anything
        through it, which may result in various exceptions being thrown. Therefore we
        try and parse as much as possible with argparse to make sure the inputs are
        sensible.

        Returns (variable, nx_expr, ny_expr, degree); nx/ny are left as strings to be
        evaluated in invoke (in the current frame).
        """
        args = DisplayCGVector._make_argument_parser().parse_args(
            gdb.string_to_argv(arg)
        )
        return args.variable, args.nx, args.ny, args.degree

    def invoke(self, arg, from_tty):
        """Invoke the command."""
        # Parse the command line arguments
        variable, nx_expr, ny_expr, degree = self._parse_arguments(arg)

        # Get the variable into a gdb.Value
        val = gdb.parse_and_eval(variable)

        # CGVector is a flat Eigen column vector;
        # In GDB `size() method tends to run into an infinite loop, so we use rows() instead.`
        n_nodes = int(gdb.parse_and_eval(f"{variable}.rows()"))

        # nx/ny are GDB expressions supplied by the user
        nx = int(gdb.parse_and_eval(nx_expr))
        ny = int(gdb.parse_and_eval(ny_expr))
        if degree < 1:
            msg = f"CG degree must be >= 1 (got {degree})"
            raise ValueError(msg)

        n_nodes_x = degree * nx + 1
        n_nodes_y = degree * ny + 1
        expected = n_nodes_x * n_nodes_y
        if expected != n_nodes:
            msg = (
                f"(degree*nx+1)*(degree*ny+1) = ({degree}*{nx}+1)*({degree}*{ny}+1) "
                f"= {expected} does not match the CGVector size ({n_nodes})"
            )
            raise ValueError(msg)

        # Get the info about the underlying Eigen buffer (CGVector is a column vector)
        storage_info = val["m_storage"]
        ptr = storage_info["m_data"]

        # Establish a numpy array on the data (flat, single component)
        data = np.array(read_to_list(ptr, n_nodes, 1))

        cg_field_plot(data, nx, ny, degree)


DisplayMArray()
DisplayDGVector()
DisplayCGVector()
