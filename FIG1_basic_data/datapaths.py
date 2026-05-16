from pathlib import Path

NCBS_locations = {
                        "figure_raw_material_location" : Path(r"../../DATA_Nov2025"),
                        "paper_figure_export_location" : Path(r"."),
                        "data_path_FS"                 : Path(r"../../DATA_Nov2025"),
                        "data_path_LTM"                : Path(r"../../DATA_Nov2025"),
                        "data_path_grid"               : Path(r"../../DATA_Nov2025"),
                        "data_path_analysed"           : Path(r"../../DATA_Nov2025"),
                        "project_path_root"            : Path(r".")
                    }

onedrive_locations = {
                        "figure_raw_material_location" : Path(r"C:\\Users\\Aditya\\OneDrive\\NCBS\\Lab\Projects\EI_Dynamics\Analysis\paper_figures\submission"),
                        "paper_figure_export_location" : Path(r"C:\\Users\\Aditya\\OneDrive\\NCBS\\Lab\Projects\EI_Dynamics\Analysis\paper_figures\submission"),
                        "data_path_FS"                 : Path(r"C:\\Users\\Aditya\\OneDrive\\NCBS\\Lab\Projects\EI_Dynamics\Analysis\parsed_data\FreqSweep"),
                        "data_path_LTM"                : Path(r"C:\\Users\\Aditya\\OneDrive\\NCBS\\Lab\Projects\EI_Dynamics\Analysis\parsed_data\LTMRand"),
                        "data_path_grid"               : Path(r"C:\\Users\\Aditya\\OneDrive\\NCBS\\Lab\Projects\EI_Dynamics\Analysis\parsed_data\Grid"),
                        "data_path_analysed"           : Path(r"C:\\Users\\Aditya\\OneDrive\\NCBS\\Lab\Projects\EI_Dynamics\Analysis\parsed_data\second_order"),
                        "project_path_root"            : Path(r"C:\\Users\\Aditya\\OneDrive\\NCBS")
                    }

# check if NCBS storage is available
if NCBS_locations["project_path_root"].exists():
    print("working on NCBS storage drive")
    location = NCBS_locations
else:
    print("working on local onedrive")
    location = onedrive_locations

