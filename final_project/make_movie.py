import imageio.v2 as imageio
import glob
import os

def stitch_frames_to_mp4(output_filename='heat_diffusion.mp4', fps=15):
    # Gather all PNG files matching the naming convention and sort them sequentially
    file_list = sorted(glob.glob("frame_*.png"))
    
    if not file_list:
        print("No frames found.")
        return

    print(f"Found {len(file_list)} frames. Stitching...")

    # Open the video writer
    writer = imageio.get_writer(output_filename, fps=fps, macro_block_size=None)
    
    for filename in file_list:
        image = imageio.imread(filename)
        writer.append_data(image)
        
    writer.close()
    print(f"Video saved as {output_filename}")

if __name__ == "__main__":
    stitch_frames_to_mp4(fps=30)