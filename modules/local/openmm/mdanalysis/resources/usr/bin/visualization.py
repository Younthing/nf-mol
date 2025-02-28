#!/usr/bin/env python3
from pathlib import Path
from typing import List, Optional

import matplotlib.pyplot as plt
import MDAnalysis as mda
import nglview as nv
import numpy as np
import pandas as pd
from MDAnalysis.analysis import align, diffusionmap, rms


class MDAnalysisBase:
    """Base class for molecular dynamics analysis."""

    def __init__(self, topology_file: str, trajectory_file: str):
        """
        Initialize molecular dynamics analysis base class.

        Args:
            topology_file: Path to topology file (.pdb format)
            trajectory_file: Path to trajectory file (.xtc format)

        Raises:
            FileNotFoundError: If topology or trajectory file does not exist
            RuntimeError: If creating the Universe object fails
        """
        if not Path(topology_file).exists() or not Path(trajectory_file).exists():
            raise FileNotFoundError("Topology file or trajectory file does not exist")

        try:
            self.universe = mda.Universe(topology_file, trajectory_file)
        except Exception as e:
            raise RuntimeError(f"Failed to create Universe object: {str(e)}")


class MDVisualizer(MDAnalysisBase):
    """Class for visualizing molecular dynamics trajectories."""

    def align_and_visualize(self, selection: str = "protein") -> nv.NGLWidget:
        """
        Align molecular dynamics trajectory and return visualization view.

        Args:
            selection: Atom selection string for alignment (default: "protein")

        Returns:
            nv.NGLWidget: NGL view object
        """
        # Reset to first frame
        self.universe.trajectory[0]

        # Align trajectory
        alignment = align.AlignTraj(
            mobile=self.universe,
            reference=self.universe,
            select=selection,
            in_memory=True,
        )
        alignment.run()

        # Return view object
        view = nv.show_mdanalysis(self.universe)
        return view


class RMSDAnalyzer(MDAnalysisBase):
    """Class for RMSD analysis of molecular dynamics trajectories."""

    def calculate_rmsd(
        self, selection1: str, selection2: Optional[List[str]] = None
    ) -> pd.DataFrame:
        """
        Calculate RMSD for selected atom groups.

        Args:
            selection1: Primary atom group selection string, also used for alignment
            selection2: List of additional atom group selection strings

        Returns:
            pandas.DataFrame: DataFrame containing RMSD data over time for selected atom groups
        """
        # Reset to first frame
        self.universe.trajectory[0]

        # Run RMSD analysis
        rmsd_analysis = rms.RMSD(
            self.universe, self.universe, select=selection1, groupselections=selection2
        )
        rmsd_analysis.run()

        # Create DataFrame from results
        columns = [selection1, *(selection2 or [])]
        rmsd_df = pd.DataFrame(
            np.round(rmsd_analysis.results.rmsd[:, 2:], 2), columns=columns
        )
        rmsd_df.index.name = "frame"
        return rmsd_df

    def plot_rmsd(
        self, rmsd_df: pd.DataFrame, title: str = "RMSD of protein and ligand"
    ) -> plt.Axes:
        """
        Create RMSD plot.

        Args:
            rmsd_df: RMSD data DataFrame
            title: Plot title

        Returns:
            plt.Axes: matplotlib plot object
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        rmsd_df.plot(ax=ax)
        ax.set_title(title)
        ax.set_xlabel("Frame")
        ax.set_ylabel("RMSD (Å)")
        ax.legend(loc="best")
        plt.tight_layout()
        return ax

    def calculate_frame_rmsd(self, selection: str) -> np.ndarray:
        """
        Calculate RMSD matrix between all frames.

        Args:
            selection: Atom group selection string

        Returns:
            np.ndarray: RMSD distance matrix
        """
        pairwise_rmsd = diffusionmap.DistanceMatrix(self.universe, select=selection)
        pairwise_rmsd.run()
        return pairwise_rmsd.results.dist_matrix

    def plot_frame_rmsd(self, protein_selection: str, ligand_name: str) -> plt.Figure:
        """
        Plot heatmap of frame-to-frame RMSD for protein and ligand.

        Args:
            protein_selection: Protein selection string
            ligand_name: Ligand residue name

        Returns:
            plt.Figure: matplotlib figure object
        """
        dist_matrix_protein = self.calculate_frame_rmsd(protein_selection)
        dist_matrix_ligand = self.calculate_frame_rmsd(f"resname {ligand_name}")

        max_dist = max(np.amax(dist_matrix_ligand), np.amax(dist_matrix_protein))

        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle("Frame-to-Frame RMSD Analysis")

        img1 = ax[0].imshow(dist_matrix_protein, cmap="viridis", vmin=0, vmax=max_dist)
        ax[0].set_title("Protein")
        ax[0].set_xlabel("Frame")
        ax[0].set_ylabel("Frame")

        img2 = ax[1].imshow(dist_matrix_ligand, cmap="viridis", vmin=0, vmax=max_dist)
        ax[1].set_title("Ligand")
        ax[1].set_xlabel("Frame")

        cbar = fig.colorbar(
            img1, ax=ax, orientation="horizontal", fraction=0.1, label="RMSD (Å)"
        )
        plt.tight_layout()
        return fig


class PocketVisualizer(MDAnalysisBase):
    """Class for visualizing binding pockets in molecular dynamics trajectories."""

    def visualize_binding_pocket(
        self, ligand_name: str, distance: float = 5.0, representation: str = "licorice"
    ) -> nv.NGLWidget:
        """
        Visualize ligand binding pocket.

        Args:
            ligand_name: Ligand residue name
            distance: Radius around ligand (Angstroms)
            representation: Visualization representation style

        Returns:
            nv.NGLWidget: NGL view object
        """
        # Find residues in the binding pocket
        pocket_atoms = self.universe.select_atoms(
            f"(around {distance} resname {ligand_name}) and protein"
        )
        pocket_resids = set(pocket_atoms.resids)

        # Create view
        view = nv.show_mdanalysis(self.universe)

        # Add representations
        view.add_representation(
            representation,
            selection=f"protein and resid {' '.join(str(x) for x in pocket_resids)}",
        )
        view.add_representation("ball+stick", selection=f"resname {ligand_name}")

        # Center view on ligand
        view.center(f"resname {ligand_name}")
        return view


def run_analysis(
    topology: str,
    trajectory: str,
    ligand_name: str = "LIG",
    output_dir: Optional[str] = None,
) -> dict:
    """
    Run a complete MD analysis workflow.

    Args:
        topology: Path to topology file
        trajectory: Path to trajectory file
        ligand_name: Ligand residue name
        output_dir: Directory to save output files (default: current directory)

    Returns:
        dict: Dictionary containing analysis results and output paths
    """
    # Set up output directory
    output_path = Path(output_dir) if output_dir else Path.cwd()
    output_path.mkdir(parents=True, exist_ok=True)

    results = {}

    try:
        # Visualization analysis
        visualizer = MDVisualizer(topology, trajectory)
        view = visualizer.align_and_visualize()
        view.render_image(trim=True, factor=2, transparent=True)
        # Save visualization
        view_path = output_path / "aligned_md.html"
        nv.write_html(str(view_path), [view], frame_range=[0, view.frame - 1])
        results["aligned_trajectory_path"] = view_path

        # RMSD analysis
        analyzer = RMSDAnalyzer(topology, trajectory)
        rmsd_df = analyzer.calculate_rmsd(
            "backbone", ["protein", f"resname {ligand_name}"]
        )
        # Save RMSD data
        rmsd_csv_path = output_path / "rmsd_results.csv"
        rmsd_df.to_csv(rmsd_csv_path)
        results["rmsd_data_path"] = rmsd_csv_path
        # Plot and save RMSD chart
        ax = analyzer.plot_rmsd(rmsd_df)
        rmsd_plot_path = output_path / "rmsd_plot.png"
        plt.savefig(rmsd_plot_path, dpi=300)
        plt.close()
        results["rmsd_plot_path"] = rmsd_plot_path

        # Calculate and plot frame-to-frame RMSD heatmap
        frame_rmsd_fig = analyzer.plot_frame_rmsd("protein", ligand_name)
        frame_rmsd_path = output_path / "frame_rmsd_plot.png"
        frame_rmsd_fig.tight_layout()
        frame_rmsd_fig.savefig(frame_rmsd_path, dpi=300)
        plt.close()
        results["frame_rmsd_path"] = frame_rmsd_path

        # Binding pocket visualization
        # Binding pocket visualization
        pocket_viewer = PocketVisualizer(topology, trajectory)
        pocket_view = pocket_viewer.visualize_binding_pocket(ligand_name)
        pocket_view_path = output_path / "pocket_view.html"
        pocket_view.center(selection=f"resname {ligand_name}")
        nv.write_html(
            str(pocket_view_path), [pocket_view], frame_range=[0, pocket_view.frame - 1]
        )
        results["pocket_view_path"] = pocket_view_path

        return results

    except Exception as e:
        print(f"Error during analysis: {str(e)}")
        raise


# Usage example:
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run MD analysis and visualization.")
    parser.add_argument("topology", type=str, help="Path to topology file (.pdb)")
    parser.add_argument("trajectory", type=str, help="Path to trajectory file (.xtc)")
    parser.add_argument("ligand_name", type=str, help="Ligand residue name")
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        default="output",
        help="Output directory (default: output)",
    )

    args = parser.parse_args()

    results = run_analysis(
        args.topology, args.trajectory, args.ligand_name, args.output_dir
    )

    print("Analysis complete. Output files:")
    for key, path in results.items():
        print(f"- {key}: {path}")
