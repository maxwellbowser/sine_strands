"""Written using Python 3.11"""

import tkinter as tk
from tkinter import ttk
from tkinter import filedialog as fd
from tkinter.messagebox import showinfo

import sys
import os
from pathlib import Path
import time
import pandas as pd
import numpy
from scipy import optimize
import matplotlib.pyplot as plt
import re
import json
from scipy.ndimage import gaussian_filter

pd.options.mode.chained_assignment = None
# This script will take DAT Files as input, separate out the frequencies programmed and
# Collect/average the waves


# Although there may be some slip-ups, Actual refers to what the sinusoidal machine is doing
# and Calc refers to the calcuated values for waves


def three_sin_calc(Hz_list, sampling_rate):
    actual_timing = []
    # for lower sampling rates, the timing interval is larger
    # Usually sampling rate is 10khz to 5khz
    # so timing_interval is ~ .1-.2
    timing_interval = 10_000 / sampling_rate

    for frequency in Hz_list:
        # this is the period of one wave in ms
        wave_timing = 1_000 / float(frequency)

        # Making sure that numbers are rounded to the tenth and everything is always rounded down!

        """
        if round(wave_timing, 1) > wave_timing:
            wave_timing = round(wave_timing, 1) - (timing_interval)

        else:
            wave_timing = round(wave_timing, 1)

        """
        wave_timing = round(wave_timing, 1)

        # This is for 5 kHz compatiability,
        # to make sure that the cut lines for sine waves lands on even time points
        # (0.2, 0.4, 0.6...ect)
        if (round(wave_timing, 1) * 10) % timing_interval != 0:
            wave_timing = round(wave_timing - 0.1, 1)

        triple_wave_timing = [float(frequency), round(wave_timing, 1) * 3]
        actual_timing.append(triple_wave_timing)

    return actual_timing


def fit_sin(tt, yy):
    """
    Fit sin to the input time sequence, and return fitting parameters:
    "amp",
    "omega",
    "phase",
    "offset",
    "freq",
    "period", and
    "fitfunc"
    Source: https://stackoverflow.com/questions/16716302/how-do-i-fit-a-sine-curve-to-my-data-with-pylab-and-numpy, Author: unsym
    """
    tt = numpy.array(tt)
    yy = numpy.array(yy)
    ff = numpy.fft.fftfreq(len(tt), (tt[1] - tt[0]))  # assume uniform spacing
    Fyy = abs(numpy.fft.fft(yy))

    # These are just taking a guess at where the algorithm should start
    # Nothing super crazy, but they use fft to pick the freq
    guess_freq = abs(
        ff[numpy.argmax(Fyy[1:]) + 1]
    )  # excluding the zero frequency "peak", which is related to offset
    guess_amp = numpy.std(yy) * 2.0**0.5
    guess_offset = numpy.mean(yy)
    guess = numpy.array([guess_amp, 2.0 * numpy.pi * guess_freq, 0.0, guess_offset])

    def sinfunc(t, A, w, p, c):
        # This is suspicious because the formula should be:
        # A * numpy.sin(w * (t + p)) + c
        # but when I do that, the curve does not fit to the data
        return A * numpy.sin(w * t + p) + c

    popt, pcov = optimize.curve_fit(sinfunc, tt, yy, p0=guess, maxfev=2000)
    A, w, p, c = popt
    f = w / (2.0 * numpy.pi)
    fitfunc = lambda t: A * numpy.sin(w * t + p) + c

    return {
        "amp": A,
        "omega": w,
        "phase": p,
        "offset": c,
        "freq": f,
        "period": 1.0 / f,
        "fitfunc": fitfunc,
        "maxcov": numpy.max(pcov),
        "rawres": (guess, popt, pcov),
    }


def main():
    """
    Sine Strands is a program for automating the spliting
    and averaging of DAT files, outputing Elastic and Viscous Moduli

    """

    def select_files():
        filetypes = (("All Files", "*.*"),)

        global filepath

        filepath = fd.askopenfilenames(
            title="Open files", initialdir=r"C:\Users\Desktop", filetypes=filetypes
        )

        root.destroy()

    filtering_amount = 25

    # The values of pCa for Relax and Active were given to me by Gerrie, may change in the future
    # Just directly from our Excel sheet
    pCa_dict = {
        "100": 4.333,
        "99": 4.563,
        "98": 4.750,
        "97": 4.894,
        "96": 5.008,
        "95": 5.102,
        "94": 5.182,
        "93": 5.250,
        "92": 5.311,
        "91": 5.366,
        "90": 5.416,
        "89": 5.462,
        "88": 5.504,
        "87": 5.543,
        "86": 5.580,
        "85": 5.615,
        "84": 5.648,
        "83": 5.679,
        "82": 5.709,
        "81": 5.738,
        "80": 5.765,
        "79": 5.792,
        "78": 5.818,
        "77": 5.843,
        "76": 5.867,
        "75": 5.890,
        "74": 5.913,
        "73": 5.935,
        "72": 5.957,
        "71": 5.978,
        "70": 5.999,
        "69": 6.020,
        "68": 6.040,
        "67": 6.060,
        "66": 6.079,
        "65": 6.098,
        "64": 6.117,
        "63": 6.136,
        "62": 6.154,
        "61": 6.173,
        "60": 6.191,
        "59": 6.209,
        "57": 6.227,
        "56": 6.262,
        "55": 6.280,
        "54": 6.297,
        "53": 6.315,
        "52": 6.332,
        "51": 6.350,
        "50": 6.367,
        "49": 6.384,
        "48": 6.402,
        "47": 6.419,
        "46": 6.437,
        "45": 6.454,
        "44": 6.472,
        "43": 6.489,
        "42": 6.507,
        "41": 6.525,
        "40": 6.543,
        "39": 6.561,
        "38": 6.579,
        "37": 6.598,
        "36": 6.617,
        "35": 6.636,
        "34": 6.655,
        "33": 6.674,
        "32": 6.694,
        "31": 6.714,
        "30": 6.735,
        "28": 6.777,
        "27": 6.799,
        "26": 6.821,
        "25": 6.844,
        "24": 6.868,
        "23": 6.892,
        "22": 6.917,
        "21": 6.942,
        "20": 6.969,
        "19": 6.997,
        "18": 7.025,
        "17": 7.056,
        "16": 7.087,
        "15": 7.120,
        "14": 7.155,
        "13": 7.193,
        "12": 7.232,
        "11": 7.275,
        "10": 7.321,
    }

    # Gui Creation:

    # Create the main window
    root = tk.Tk()
    root.title("Sine Strands")

    # Create and configure the frame
    frame = ttk.Frame(root, padding=10)
    frame.grid(column=0, row=0)

    # Making variables
    filtering_var_Fin = tk.BooleanVar(value=False)
    filtering_var_Lin = tk.BooleanVar(value=False)
    plots_var = tk.BooleanVar(value=False)
    folder_name = tk.StringVar()
    filtering_strength = tk.IntVar(value=filtering_amount)

    # Create and place the label for folder name
    folder_name_label = ttk.Label(frame, text="Folder Name:")
    folder_name_label.grid(column=0, row=0, sticky="w", padx=10, pady=10)

    # Create and place the entry field for folder name
    folder_name_entry = tk.Entry(frame, textvariable=folder_name)
    folder_name_entry.grid(column=1, row=0, padx=10, pady=10)

    # Create and place the checkbox for enabling filtering
    filtering_Fin_checkbox = ttk.Checkbutton(
        frame, text="Filtering (Fin)", variable=filtering_var_Fin
    )
    filtering_Fin_checkbox.grid(column=0, row=2, sticky="w", padx=10, pady=10)

    filtering_Lin_checkbox = ttk.Checkbutton(
        frame, text="Filtering (Fin)", variable=filtering_var_Lin
    )
    filtering_Lin_checkbox.grid(column=0, row=3, sticky="w", padx=10, pady=10)

    plots_checkbox = ttk.Checkbutton(frame, text="Em vs. Vm Plots", variable=plots_var)
    plots_checkbox.grid(column=1, row=2, sticky="w", padx=10, pady=10)

    # Create and place the label and entry for filter strength
    filter_strength_label = ttk.Label(frame, text="Filter Strength:")
    filter_strength_label.grid(column=0, row=1, sticky="w", padx=10, pady=10)

    filter_strength_entry = tk.Entry(frame, textvariable=filtering_strength)
    filter_strength_entry.grid(column=1, row=1, padx=10, pady=10)

    # Create and place the button to open the files
    create_button = ttk.Button(frame, text="Open Files", command=select_files)
    create_button.grid(column=1, row=3, columnspan=2, padx=10, pady=10)

    # Run the tkinter main loop
    root.mainloop()

    # Making the first folder for all of the files
    current_dir = Path(os.getcwd())
    new_dir = current_dir / folder_name.get()

    try:
        os.mkdir(new_dir)
    except FileExistsError:
        print(
            "Sorry, it looks like that folder already exists! Please delete it or try a different filename"
        )
        time.sleep(2)
        main()

    os.chdir(new_dir)

    # This is used to keep track if any files didn't have a cell area
    no_cross_area = False
    cross_errors = {"Non-Scaled Files:": []}

    # Keeping track of cell forces, for each file
    force_dict = {"Avg. Basal Force": [], "Filename": []}

    for file in filepath:
        filename = os.path.basename(file)

        # Folder Creation for that file
        current_dir = Path(os.getcwd())
        new_dir = current_dir / filename[:-4]
        os.mkdir(new_dir)
        os.chdir(new_dir)

        # Using Regex for finding SL
        # just looks for 3 digits in a string
        pattern = r"\d{3}"
        matches = re.findall(pattern, filename[-7:-4])

        # If there are 3 digits in a row at the end of the string, make that the SL
        if len(matches) == 1:
            SL = matches[0]
            SL = int(SL)

        # Otherwise, make the SL 225 (default)
        else:
            SL = 225

        # Looking for the conditions in the filename
        if "treated" in filename.lower():
            condition = "Treated"

        elif "untreated" in filename.lower():
            condition = "Untreated"

        else:
            condition = "None"

        # If the file name has "relax", then there is 0 % activating
        if "relax" in filename.lower():
            pCa = 9.0

        # Else, if there is "100" in the filename, set 100 %
        # with Regex, the other option was to look for 2 or 3 digits in the filename,
        # but this caused names if there was also a different SL #
        elif "100" in filename:
            pCa = pCa_dict["100"]

        else:
            # This gives all 2 digit numbers, preceeded by the letter "d" (treateD or untreateD)
            # so it doesn't capture any 3 digit codes like if SL is present at the end of the filename
            percent_activating = re.findall(r"d(\d{2})", filename)

            # Find the pCa value matching to the % of activating (gotten from Gerries Excel calculator)
            pCa = pCa_dict[percent_activating[0]]

        # Getting Mouse ID
        # Regex from Claude.ai
        # The numbers after wt or ko and before un (untreated) or tr (treated)
        # Here's a breakdown of the pattern:
        #     (?:ko|wt): This non-capturing group matches either "ko" or "wt".
        #     (\d+): This matches and captures the one or more digits (numbers).
        #     (?:un|tr): This non-capturing group matches either "un" or "tr".
        mouse_ID_pattern = r"(?:ko|wt)(\d+)(?:un|tr)"
        matches = re.findall(mouse_ID_pattern, filename.lower())

        # First match should be the mouse ID (should only be one result)
        mouse_ID = matches[0]

        # TBH I got this from stackoverflow
        try:
            DAT_df = pd.DataFrame([i.strip().split() for i in open(file).readlines()])

        # Check if file is corrupted, if it is then skip it!
        # Also give some time for the person to read the filename
        except UnicodeDecodeError:
            print(f"Uh Oh, Corrupted File! Please check: {filename}")
            time.sleep(2)
            continue

        # This is NOT adaptive, need to fix in future
        sampling_rate = float(DAT_df.iat[2, 3])
        sampling_scale = 10_000 / sampling_rate

        # For the files that Gerrie creates, he stores the fiber thickness and width in the comments, which is in the 0
        # column and 4/5th position rows. This code below calculates the cross sectional area of the fiber,
        # or if those numbers are not present, then it sets the area to None
        try:
            fiber_info = (float(DAT_df.iat[4, 0]), float(DAT_df.iat[5, 0]))
            cross_area = (fiber_info[0] * 0.5) * (fiber_info[1] * 0.5) * numpy.pi

        except ValueError:
            # This is if Gerrie added "width" and "length" as spacers, if they're not there then idk
            try:
                fiber_info = (float(DAT_df.iat[4, 1]), float(DAT_df.iat[5, 1]))
                cross_area = (fiber_info[0] * 0.5) * (fiber_info[1] * 0.5) * numpy.pi

            except ValueError:
                cross_area = None

        # DAT_df.to_csv(f"{os.path.basename(filename)}.csv")

        # This is finding where the parameters start and end
        # Looks for "time", where the second and third instance are the boundries for the timings
        parameter_locs = DAT_df.loc[DAT_df[0] == "Time"].index.tolist()
        parameter_locs[2] -= 2

        # Using parameter_locs to actually separate the values from the full df
        DAT_parameters = DAT_df.iloc[parameter_locs[1] : parameter_locs[2], :8]
        DAT_parameters.reset_index(drop=True, inplace=True)

        # Next, we need to find out if the timings are ms or s
        # I'm just looking in column 7 for s/ms by looking halfway down the list
        # I'm going to convert everything to ms, so if its in seconds, we set it to 1000ms
        if DAT_parameters.iat[len(DAT_parameters) // 2, 7] == "s":
            s_vs_ms = 1000

        elif DAT_parameters.iat[len(DAT_parameters) // 2, 7] == "ms":
            s_vs_ms = 1

        # If nothing is found, then we just look one more down, incase there's spacing
        else:
            if DAT_parameters.iat[(len(DAT_parameters) // 2) + 1, 7] == "s":
                s_vs_ms = 1000

            elif DAT_parameters.iat[(len(DAT_parameters) // 2) + 1, 7] == "ms":
                s_vs_ms = 1

        # Getting the Hz and length of Hz application times
        hz_plus_timings = DAT_parameters.iloc[:, [2, 6]]

        # Getting the indexes of where the Hz are (basically every other row)
        locations_of_hz = DAT_parameters.loc[DAT_parameters[3] == "Hz"].index.tolist()

        # separating out the rows which only have "Hz" in the 3rd column
        # basically skipping all step instructions
        prgrm_timing = hz_plus_timings.iloc[locations_of_hz]

        # Changing labels to more descriptive ones
        prgrm_timing.columns = ["Hz", "Time"]

        # Scaling the timings to ms, if they are set as sec
        # (basically multiply by 1000 (timing in s) or 1 (timing already in ms))
        prgrm_timing.loc[:, "Time"] = (
            prgrm_timing.loc[:, "Time"].astype("float").mul(s_vs_ms).round(1)
        )

        prgrm_timing.reset_index(drop=True, inplace=True)

        # self-explanitory, getting all of the frequencies as a list
        frequencies = list(prgrm_timing.loc[:, "Hz"])

        # Calculating the perfect (scaled down) frequencies
        trimmed_timings = three_sin_calc(frequencies, sampling_rate)

        calc_timings_df = pd.DataFrame(trimmed_timings, columns=["Hz", "Calc"])

        actual_times_df = pd.DataFrame(prgrm_timing.values, columns=["Hz", "Actual"])
        actual_times_df["Actual"] = pd.to_numeric(actual_times_df["Actual"])

        # Now we have the calculated times & machine command times for 3 sine waves at a given frequency
        # order is {0:"Hz", 1:"Calc", 2:"Actual"}
        combined_timing_df = pd.concat(
            [calc_timings_df, actual_times_df.iloc[:, 1]], axis="columns"
        )

        # Preparing the data to be analyzed
        column_names = [
            "Time",
            "Lin",
            "Lout",
            "Fin",
            "Fout",
            "SL",
        ]

        fiber_data = pd.DataFrame(
            DAT_df.iloc[(parameter_locs[2] + 3) :, [0, 1, 2, 3, 4, 7]]
        )
        fiber_data.reset_index(drop=True, inplace=True)
        # fiber_data.drop(fiber_data.columns[9:], axis=1, inplace=True)
        fiber_data.columns = column_names

        outputs = {
            "Freq (Hz)": [],
            "Fin_Amplitude": [],
            "Lin_Amplitude": [],
            "Fin_PhaseShift": [],
            "Lin_PhaseShift": [],
            "CrossArea": [],
            "Em (kPa)": [],
            "Vm (kPa)": [],
        }

        # At this point in the script, the data is separated, indexes are reset,
        # the column names have been changed, and we are ready to do some math (yaaaayyy)
        # print(fiber_data)
        # print(combined_timing_df)
        # print(fiber_data.tail())

        # Setting the time point where the first sine wave begins
        finding_start = DAT_parameters.loc[
            DAT_parameters[1] == "Length-Sine"
        ].index.to_list()

        # First instance of "Length-Sine" is the line where the data start
        start = DAT_parameters.iloc[finding_start[0], 0]

        scaling_factor = 10 / sampling_scale

        start = int(
            float(start) * scaling_factor
        )  # Since the setup takes 200ms to start going
        spacing = 10 * scaling_factor

        # This is the average first 200 values of the force, before the sinsusoidal perturbation starts...
        first_force = (
            fiber_data["Fin"].iloc[int(start) : int(start) + 200].astype("float").mean()
        )
        force_dict["Avg. Basal Force"].append(first_force)
        force_dict["Filename"].append(filename[:-4])

        for freq in range(0, len(combined_timing_df)):
            end = int(
                start + (combined_timing_df["Actual"].iloc[freq] * scaling_factor)
            )

            full_data_df = fiber_data.iloc[int(start) : int(end)]
            full_data_df.reset_index(drop=True, inplace=True)
            full_data_df = full_data_df.astype(
                {
                    "Time": "float",
                    "Lin": "float",
                    "Lout": "float",
                    "Fin": "float",
                    "Fout": "float",
                    "SL": "float",
                }
            )

            # Subtracting the time per frequency of actual - calculated
            # If the timing for a frequency divides evenly into thirds,
            # timing difference will be zero. If it doesn't, then timing_difference gives the
            timing_difference = scaling_factor * round(
                combined_timing_df.iat[freq, 2] - combined_timing_df.iat[freq, 1], 1
            )

            # This is a weird thing but ya know if it works, it works
            # basically if there is a timing difference and if there is a remainder when
            # dividing
            if (
                timing_difference != 0
                and sampling_rate % combined_timing_df.iat[freq, 0] != 0
            ):
                full_data_df = full_data_df.iloc[: -int(timing_difference)]

            # This here is from stack overflow, creating an array of numbers from 0 to the time of frequency, that has a length of the Timepoints.
            # Basically making an array to copy the # of ms passed on the file
            x_data = numpy.linspace(
                0, combined_timing_df.iat[freq, 2], len(full_data_df.Time)
            )

            # x_data = full_data_df.Time
            # base = peakutils.baseline(y_data)
            # y_data = y_data - base
            if filtering_var_Fin.get():
                # Filtering the Force
                full_data_df["Fin"] = gaussian_filter(
                    full_data_df["Fin"], filtering_strength.get()
                )
            elif filtering_var_Lin.get():
                # Filtering the Force
                full_data_df["Lin"] = gaussian_filter(
                    full_data_df["Lin"], filtering_strength.get()
                )

            # Fitting sine curves
            fin_res = fit_sin(x_data, full_data_df.Fin)
            lin_res = fit_sin(x_data, full_data_df.Lin)

            # Code for saving plots of every freqency and the overlaid Sine fit
            # adds a lot of time and slows down everything
            """
            plt.plot(x_data, full_data_df.Fin)
            plt.plot(
                x_data,
                fin_res["fitfunc"](x_data),
                "r-",
                label="y fit curve",
                linewidth=2,
            )

            plt.savefig(f"fin_{combined_timing_df.iloc[freq, 0]}.png")
            plt.close()

            plt.plot(x_data, full_data_df.Lin)
            plt.plot(
                x_data,
                lin_res["fitfunc"](x_data),
                "r-",
                label="y fit curve",
                linewidth=2,
            )

            plt.savefig(f"lin_{combined_timing_df.iloc[freq, 0]}.png")
            plt.close()
            """

            # Try to scale by the cell's cross-sectional area, unless there is none found, then return the raw amplitude of force
            try:
                fin_amp = fin_res["amp"] / cross_area

            except TypeError:
                # Need to add a thing that keeps track of non-scaled values
                no_cross_area = True
                fin_amp = fin_res["amp"]
                cross_errors["Non-Scaled Files:"].append(filename)

            # full_data_df.to_csv(f"Output_{combined_timing_df.iloc[freq, 0]}_Hz.csv")
            norm_phase = fin_res["phase"] - lin_res["phase"]

            elastic_calc = (fin_amp / lin_res["amp"]) * numpy.cos(norm_phase)
            viscous_calc = (fin_amp / lin_res["amp"]) * numpy.sin(norm_phase)

            outputs["Freq (Hz)"].append(combined_timing_df.iloc[freq, 0])
            outputs["Fin_Amplitude"].append(fin_amp)
            outputs["Lin_Amplitude"].append(lin_res["amp"])
            outputs["Fin_PhaseShift"].append(fin_res["phase"])
            outputs["Lin_PhaseShift"].append(lin_res["phase"])
            outputs["CrossArea"].append(cross_area)
            outputs["Em (kPa)"].append(elastic_calc)
            outputs["Vm (kPa)"].append(viscous_calc)

            # Remember! micro (u) sec not milli (m) sec
            # print(start, end)
            start = end + spacing

        # os.chdir("../")
        summary_df = pd.DataFrame.from_dict(outputs)

        if "concat_df" in locals():
            # 1 mN/mm2 = 1 kPa
            adding_df = summary_df.loc[:, ["Freq (Hz)", "Vm (kPa)", "Em (kPa)"]]

            # Adding Filename
            adding_df.insert(
                0,
                "Filename",
                filename[:-4],
                allow_duplicates=True,
            )

            # Adding HeartSample (for MATLAB)
            adding_df.insert(
                0,
                "HeartSample",
                "Mouse" + mouse_ID,
                allow_duplicates=True,
            )

            # Adding SL (for MATLAB)
            adding_df.insert(
                0,
                "SL",
                SL,
                allow_duplicates=True,
            )

            # Adding Treatment (for MATLAB)
            #
            adding_df.insert(
                0,
                "Treatment",
                condition,
                allow_duplicates=True,
            )

            # Adding pCa (for MATLAB)
            # !! Skeletal are all 76 !!
            # pCa, relaxed
            # 76 = 5.87 76% activating free ca is 1.3575
            # 78 = 5.82 78% free Ca is 1.512
            # relaxed long is nothing
            # relaxed short is the useful one
            # Relaxed = pCa 9
            adding_df.insert(
                0,
                "pCa",
                pCa,
                allow_duplicates=True,
            )

            # Adding HashCode (for MATLAB) Used "sinusoid" for the CRC32 HashCode 9efeb779 (random)
            adding_df.insert(
                0,
                "HashCode",
                "9efeb779",
                allow_duplicates=True,
            )

            concat_df = pd.concat([concat_df, adding_df])

        else:
            # Making first set of data for the full data
            concat_df = summary_df.loc[:, ["Freq (Hz)", "Vm (kPa)", "Em (kPa)"]]
            # Adding Filename
            concat_df.insert(
                0,
                "Filename",
                filename[:-4],
                allow_duplicates=True,
            )

            # Adding HeartSample (for MATLAB)
            concat_df.insert(
                0,
                "HeartSample",
                "Mouse" + mouse_ID,
                allow_duplicates=True,
            )

            # Adding SL (for MATLAB)
            # This is in the file names 200/225
            concat_df.insert(
                0,
                "SL",
                SL,
                allow_duplicates=True,
            )

            # Adding Treatment (for MATLAB)
            concat_df.insert(
                0,
                "Treatment",
                condition,
                allow_duplicates=True,
            )

            # Adding pCa (for MATLAB)
            concat_df.insert(
                0,
                "pCa",
                pCa,
                allow_duplicates=True,
            )

            # Adding HashCode (for MATLAB) Used "sinusoid" for the CRC32 HashCode 9efeb779 (random)
            concat_df.insert(
                0,
                "HashCode",
                "9efeb779",
                allow_duplicates=True,
            )

        if plots_var.get():

            # Plots for Elastic vs Viscous moduli
            plt.scatter(
                summary_df["Freq (Hz)"], summary_df["Em (kPa)"], label="Elastic"
            )
            plt.scatter(
                summary_df["Freq (Hz)"], summary_df["Vm (kPa)"], label="Viscous"
            )
            plt.title("Elastic vs. Viscous Moduli")

            plt.legend(loc="upper left")
            plt.savefig(f"{filename[:-4]}.png", dpi=150)
            plt.close()

        summary_df.to_csv(f"{filename[:-4]}_Summary_Output.csv")
        os.chdir("../")

    # Making a summary csv file which contains the columns Freq (Hz), Em (kPa), and Vm (kPa)
    force_df = pd.DataFrame.from_dict(force_dict)
    try:
        # Save dfs as excel file
        with pd.ExcelWriter(f"Full_Data_Summary-{folder_name.get()}.xlsx") as writer:
            # use to_excel function and specify the sheet_name and index
            # to store the dataframe in specified sheet
            concat_df.to_excel(
                writer,
                sheet_name="Moduli Summary",
                index=0,
            )
            force_df.to_excel(
                writer,
                sheet_name="Cell Forces",
                index=0,
            )

    except UnboundLocalError:
        print("No data was selected, press enter to continue or ctrl+c to exit...")
        input()
        main()
        sys.exit()

    # Saving filenames of files that didn't get scaled, if any exist
    if no_cross_area:
        with open(f"SineStrands_output.txt", "w") as f:
            f.write(json.dumps(cross_errors, indent=4))

        showinfo(
            title="Cell Area Missing",
            message=f'Some files did not have cell measurements! Please check "SineStrands_output" for more details. ',
        )

    print("All Done!")
    folder = os.getcwd()
    os.startfile(folder)
    sys.exit()


if __name__ == "__main__":
    main()
