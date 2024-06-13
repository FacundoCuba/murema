import sys
import pysam
import matplotlib.pyplot as plt
import numpy as np

def generate_dispersion_graph(bam_file, sample_name, ref_name, mean_depth_threshold):
    try:
        # Open the BAM file
        with pysam.AlignmentFile(bam_file, "rb") as samfile:
            # Retrieve the reference length from the BAM file header
            ref_length = samfile.header.get_reference_length(ref_name)

            # Initialize variables
            positions = []
            depth_per_position = []

            # Iterate through the reads and calculate depth of coverage
            for column in samfile.pileup():
                try:
                    positions.append(column.reference_pos)
                    depth = sum(1 for read in column.pileups if not read.is_del and not read.is_refskip)
                    depth_per_position.append(depth)
                except Exception as e:
                    print(f"Error processing read at position {column.reference_pos}: {e}")

            # Create the dispersion graph
            fig, ax1 = plt.subplots(figsize=(20, 10))

            # Scatter plot for depth
            ax1.set_xlim(left=0, right=ref_length)
            ax1.set_xlabel("Position")
            ax1.set_ylabel("Depth")
            colors_by_depth = ['red' if x == 0 else 'yellow' if (0 < x < mean_depth_threshold) else 'green' for x in depth_per_position]
            ax1.scatter(positions, depth_per_position, marker='|', color=colors_by_depth, label='Vertical Depth')

            # Calculate cumulative fraction for coverage
            ax2 = ax1.twinx()
            ax2.set_ylim(bottom=0, top=1.00)
            ax2.set_ylabel("Coverage")
            mask = np.array(depth_per_position) >= mean_depth_threshold
            cumulative_coverage = np.cumsum(mask) / len(positions)
            max_coverage = cumulative_coverage[-1]
            ax2.annotate(f'Max Coverage = {max_coverage:.2f}', xy=(1.02,max_coverage), xycoords='axes fraction', fontsize=10, color='blue', va='center')
            ax2.plot(positions, cumulative_coverage, linestyle='-', color='blue', label='Horizontal Coverage')
            plt.title(f"Depth Dispersion and Coverage Proportion for {sample_name} using {ref_name} as reference")

            # Add legend
            fig.legend(loc='lower center')

            # Save the graph as an image
            output_file = f"{sample_name}.{ref_name}.png"
            plt.savefig(output_file)
            plt.close(fig)

            # Print information about the saved files
            print(f"Dispersion graph saved as: {output_file}")

    except Exception as e:
        print(f"An error occurred while processing the BAM file: {e}")

# Entry point of the script
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python grapher.py <bam_file> <sample_name> <ref_name> <mean_depth_threshold>")
        sys.exit(1)
    
    bam_file = sys.argv[1]
    sample_name = sys.argv[2]
    ref_name = sys.argv[3]
    mean_depth_threshold = int(sys.argv[4])

    generate_dispersion_graph(bam_file, sample_name, ref_name, mean_depth_threshold)
