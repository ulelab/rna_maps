"""Logging setup and runtime reporting."""

import datetime
import logging
import os
import time


def setup_logging(output_path):
    """Sets up logging to file and console."""
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"execution_{timestamp}.log"
    log_path = os.path.join(output_path, log_filename)

    logger = logging.getLogger()
    logger.handlers.clear()
    logger.setLevel(logging.INFO)

    file_handler = logging.FileHandler(log_path, mode='w')
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    file_handler.flush = file_handler.stream.flush

    start_time = time.time()

    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
    logging.getLogger('matplotlib.backends').setLevel(logging.ERROR)
    logging.getLogger('fontTools').setLevel(logging.ERROR)
    logging.getLogger('fontTools.subset').setLevel(logging.ERROR)

    logger.info(f"Starting script execution at {timestamp}")

    return log_filename, start_time, logger


def log_runtime(start_time, logger):
    """Logs the script execution time."""
    end_time = time.time()
    total_seconds = int(end_time - start_time)
    minutes, seconds = divmod(total_seconds, 60)
    runtime_str = f"{minutes}m {seconds}s" if minutes > 0 else f"{seconds}s"
    logger.info(f"Script completed in {runtime_str}.")
    for handler in logger.handlers:
        handler.flush()
