from reportfunk.funks import custom_logger as custom_logger

def log_handler(msg):
    logger = custom_logger.Logger()
    return logger.log_handler