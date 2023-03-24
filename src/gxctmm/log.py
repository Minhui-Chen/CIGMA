import logging

logger = logging.getLogger('GxCTMM')
logger.setLevel(logging.INFO)

# format
log_format = "[%(asctime)s - %(levelname)s] %(message)s"
date_format = "%Y-%m-%d %H:%M:%S"
fmt = logging.Formatter(fmt=log_format, datefmt=date_format)
