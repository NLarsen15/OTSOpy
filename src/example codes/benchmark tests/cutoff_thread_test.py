from OTSO import cutoff
import itertools
import os
import time


if __name__ == '__main__':

    stations_list = ["OULU", "ROME", "CALG"]

    # Maximum value to test: number of logical processors available
    max_value = os.cpu_count() or 1

    # Test configurations where either corenum or threadnum is 1.
    #
    # This produces:
    #   (1, 1)
    #   (1, 2)
    #   (1, 3)
    #   (2, 1)
    #   (3, 1)
    combinations = [
        (corenum, threadnum)
        for corenum in range(1, max_value + 1)
        for threadnum in range(1, max_value + 1)
        if corenum == 1 or threadnum == 1
    ]

    # Store all results
    all_results = {}

    # ============================================================
    # Run tests
    # ============================================================

    for corenum, threadnum in combinations:

        test_name = f"corenum={corenum}, threadnum={threadnum}"

        print(f"\n{'=' * 60}")
        print(f"Testing {test_name}")
        print(f"{'=' * 60}")

        results = []

        # Run each combination 5 times
        for run in range(5):

            print(
                f"  Run {run + 1}/5...",
                end=" ",
                flush=True
            )

            start_time = time.perf_counter()

            cutoff_results = cutoff(
                Stations=stations_list,
                computation_params={
                    "corenum": corenum,
                    "threadnum": threadnum
                },
                datetime_params={
                    "year": 2000,
                    "month": 1,
                    "day": 1,
                    "hour": 0
                },
                integration_params={"gyropercent":15}
            )

            elapsed = time.perf_counter() - start_time

            # DataFrame containing Ru, Rc, Rl, etc.
            df = cutoff_results[0].copy()

            # Text output
            text_output = cutoff_results[1]

            results.append({
                "dataframe": df,
                "text": text_output,
                "time": elapsed
            })

            print(f"{elapsed:.3f} seconds")

        all_results[(corenum, threadnum)] = results


    # ============================================================
    # Check whether all 5 runs produced identical results
    # ============================================================

    print("\n\n")
    print("=" * 70)
    print("RESULT SUMMARY")
    print("=" * 70)

    overall_pass = True

    for (corenum, threadnum), results in all_results.items():

        # --------------------------------------------------------
        # Compare DataFrames
        # --------------------------------------------------------

        first_df = results[0]["dataframe"]

        dataframes_same = all(
            first_df.equals(result["dataframe"])
            for result in results[1:]
        )

        # --------------------------------------------------------
        # Compare text output
        # --------------------------------------------------------

        first_text = results[0]["text"]

        text_same = all(
            first_text == result["text"]
            for result in results[1:]
        )

        # --------------------------------------------------------
        # Timing statistics
        # --------------------------------------------------------

        times = [result["time"] for result in results]

        average_time = sum(times) / len(times)
        min_time = min(times)
        max_time = max(times)

        # --------------------------------------------------------
        # Print results
        # --------------------------------------------------------

        print(
            f"\ncorenum={corenum}, threadnum={threadnum}"
        )

        print(
            f"  DataFrame results identical: {dataframes_same}"
        )

        print(
            f"  Text results identical:      {text_same}"
        )

        print(
            f"  Average time:                {average_time:.3f} s"
        )

        print(
            f"  Min time:                    {min_time:.3f} s"
        )

        print(
            f"  Max time:                    {max_time:.3f} s"
        )

        if dataframes_same and text_same:
            print("  PASS: All 5 runs are identical.")
        else:
            print("  FAIL: Results differ between runs.")
            overall_pass = False


    # ============================================================
    # Print actual results
    # ============================================================

    print("\n\n")
    print("=" * 70)
    print("RESULTS")
    print("=" * 70)

    for (corenum, threadnum), results in all_results.items():

        print(
            f"\n--- corenum={corenum}, threadnum={threadnum} ---"
        )

        print(results[0]["dataframe"])


    # ============================================================
    # Overall result
    # ============================================================

    print("\n\n")
    print("=" * 70)
    print("OVERALL TEST RESULT")
    print("=" * 70)

    if overall_pass:
        print("PASS: All configurations produced identical results.")
    else:
        print("FAIL: At least one configuration produced different results.")

        # Make GitHub Actions / CI fail if results are inconsistent
        raise SystemExit(1)