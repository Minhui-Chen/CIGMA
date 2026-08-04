import os, sys, logging

# format
log_format = "[%(asctime)s - %(levelname)s] %(message)s"
date_format = "%Y-%m-%d %H:%M:%S"
fmt = logging.Formatter(fmt=log_format, datefmt=date_format)
level = logging.DEBUG if os.environ.get('CIGMA_DEBUG') else logging.INFO

logging.basicConfig(stream=sys.stdout, format=log_format, datefmt=date_format, level=level)
logger = logging.getLogger('CIGMA')

