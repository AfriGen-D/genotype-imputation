#!/usr/bin/env python3
"""
Comprehensive error handling and logging utility for ChiPImputation pipeline.
Provides consistent error reporting and recovery mechanisms across all modules.
"""

import sys
import os
import json
import logging
import traceback
from pathlib import Path
from datetime import datetime
from functools import wraps
from typing import Any, Dict, Optional, List

# Configure logging
def setup_logging(log_file: Optional[str] = None, level: str = "INFO"):
    """Setup comprehensive logging configuration."""
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    
    handlers = [logging.StreamHandler(sys.stderr)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=getattr(logging, level),
        format=log_format,
        handlers=handlers
    )
    
    return logging.getLogger('ChiPImputation')

class PipelineError(Exception):
    """Base exception for pipeline-specific errors."""
    def __init__(self, message: str, error_code: int = 1, context: Dict[str, Any] = None):
        super().__init__(message)
        self.error_code = error_code
        self.context = context or {}
        self.timestamp = datetime.now().isoformat()

class DataError(PipelineError):
    """Error related to data processing or validation."""
    pass

class FileError(PipelineError):
    """Error related to file operations."""
    pass

class PlottingError(PipelineError):
    """Error specific to plotting operations."""
    pass

class ImputationError(PipelineError):
    """Error specific to imputation processes."""
    pass

def error_handler(func):
    """Decorator for comprehensive error handling in pipeline functions."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        logger = logging.getLogger(func.__module__)
        
        try:
            return func(*args, **kwargs)
        except PipelineError as e:
            logger.error(f"Pipeline error in {func.__name__}: {str(e)}")
            logger.error(f"Error context: {e.context}")
            write_error_report(e, func.__name__)
            sys.exit(e.error_code)
        except FileNotFoundError as e:
            logger.error(f"File not found in {func.__name__}: {str(e)}")
            write_error_report(
                FileError(f"Required file not found: {str(e)}", 2),
                func.__name__
            )
            sys.exit(2)
        except ValueError as e:
            logger.error(f"Value error in {func.__name__}: {str(e)}")
            write_error_report(
                DataError(f"Invalid data format or value: {str(e)}", 3),
                func.__name__
            )
            sys.exit(3)
        except Exception as e:
            logger.error(f"Unexpected error in {func.__name__}: {str(e)}")
            logger.error(f"Traceback: {traceback.format_exc()}")
            write_error_report(
                PipelineError(f"Unexpected error: {str(e)}", 99),
                func.__name__
            )
            sys.exit(99)
    
    return wrapper

def write_error_report(error: PipelineError, function_name: str, 
                      output_dir: str = None):
    """Write detailed error report for debugging."""
    if output_dir is None:
        output_dir = os.environ.get('CHIPIMPUTATION_ERROR_DIR', '.')
    
    error_report = {
        'timestamp': error.timestamp,
        'function': function_name,
        'error_type': error.__class__.__name__,
        'message': str(error),
        'error_code': error.error_code,
        'context': error.context,
        'traceback': traceback.format_exc()
    }
    
    # Create error report filename
    timestamp_str = datetime.now().strftime("%Y%m%d_%H%M%S")
    error_file = Path(output_dir) / f"error_{function_name}_{timestamp_str}.json"
    
    try:
        with open(error_file, 'w') as f:
            json.dump(error_report, f, indent=2)
    except:
        # If we can't write the error report, at least log it
        print(f"ERROR REPORT: {json.dumps(error_report, indent=2)}", file=sys.stderr)

def validate_input_files(*files: str):
    """Validate that input files exist and are readable."""
    for file_path in files:
        if not file_path:
            continue
        path = Path(file_path)
        if not path.exists():
            raise FileError(
                f"Input file does not exist: {file_path}",
                context={'file': file_path, 'cwd': os.getcwd()}
            )
        if not path.is_file():
            raise FileError(
                f"Path is not a file: {file_path}",
                context={'file': file_path, 'type': 'directory' if path.is_dir() else 'other'}
            )
        if not os.access(path, os.R_OK):
            raise FileError(
                f"File is not readable: {file_path}",
                context={'file': file_path, 'permissions': oct(path.stat().st_mode)}
            )

def validate_output_directory(output_dir: str, create: bool = True):
    """Validate output directory and optionally create it."""
    path = Path(output_dir)
    
    if path.exists():
        if not path.is_dir():
            raise FileError(
                f"Output path exists but is not a directory: {output_dir}",
                context={'path': output_dir}
            )
        if not os.access(path, os.W_OK):
            raise FileError(
                f"Output directory is not writable: {output_dir}",
                context={'path': output_dir}
            )
    elif create:
        try:
            path.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            raise FileError(
                f"Could not create output directory: {output_dir}",
                context={'path': output_dir, 'error': str(e)}
            )
    else:
        raise FileError(
            f"Output directory does not exist: {output_dir}",
            context={'path': output_dir}
        )

def safe_division(numerator: float, denominator: float, 
                 default: float = 0.0) -> float:
    """Safely perform division with zero check."""
    if denominator == 0:
        return default
    return numerator / denominator

def safe_percentage(value: float, total: float, decimals: int = 2) -> float:
    """Calculate percentage safely with zero check."""
    if total == 0:
        return 0.0
    return round((value / total) * 100, decimals)

def handle_empty_data(data_name: str, logger: logging.Logger = None,
                     return_default: Any = None):
    """Handle cases where data is empty or missing."""
    msg = f"No data available for {data_name}"
    if logger:
        logger.warning(msg)
    else:
        print(f"WARNING: {msg}", file=sys.stderr)
    
    return return_default

def create_fallback_plot(output_file: str, message: str, 
                        dataset: str = "Unknown", logger: logging.Logger = None):
    """Create a fallback plot when data is insufficient."""
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.text(0.5, 0.5, f"{message}\n\nDataset: {dataset}", 
            ha='center', va='center', fontsize=14,
            transform=ax.transAxes)
    ax.set_title('Data Processing Notice')
    ax.axis('off')
    
    with PdfPages(output_file) as pdf:
        pdf.savefig(fig, bbox_inches='tight')
    plt.close()
    
    if logger:
        logger.info(f"Created fallback plot: {output_file}")
    
    return output_file

def log_pipeline_stats(stats: Dict[str, Any], logger: logging.Logger = None):
    """Log pipeline statistics in a formatted way."""
    if not logger:
        logger = logging.getLogger('ChiPImputation')
    
    logger.info("=" * 60)
    logger.info("Pipeline Statistics:")
    logger.info("=" * 60)
    
    for key, value in stats.items():
        if isinstance(value, dict):
            logger.info(f"{key}:")
            for sub_key, sub_value in value.items():
                logger.info(f"  {sub_key}: {sub_value}")
        else:
            logger.info(f"{key}: {value}")
    
    logger.info("=" * 60)

def create_error_summary(error_dir: str = ".") -> Dict[str, Any]:
    """Create a summary of all errors in the pipeline run."""
    error_files = list(Path(error_dir).glob("error_*.json"))
    
    summary = {
        'total_errors': len(error_files),
        'errors_by_type': {},
        'errors_by_function': {},
        'error_details': []
    }
    
    for error_file in error_files:
        try:
            with open(error_file) as f:
                error_data = json.load(f)
                
                # Count by type
                error_type = error_data.get('error_type', 'Unknown')
                summary['errors_by_type'][error_type] = \
                    summary['errors_by_type'].get(error_type, 0) + 1
                
                # Count by function
                function = error_data.get('function', 'Unknown')
                summary['errors_by_function'][function] = \
                    summary['errors_by_function'].get(function, 0) + 1
                
                # Add to details
                summary['error_details'].append({
                    'timestamp': error_data.get('timestamp'),
                    'function': function,
                    'type': error_type,
                    'message': error_data.get('message')
                })
        except:
            continue
    
    return summary

# Context manager for temporary error handling
class ErrorContext:
    """Context manager for localized error handling."""
    
    def __init__(self, operation: str, logger: logging.Logger = None):
        self.operation = operation
        self.logger = logger or logging.getLogger('ChiPImputation')
        self.start_time = None
    
    def __enter__(self):
        self.start_time = datetime.now()
        self.logger.debug(f"Starting {self.operation}")
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        duration = (datetime.now() - self.start_time).total_seconds()
        
        if exc_type is None:
            self.logger.debug(f"Completed {self.operation} in {duration:.2f}s")
        else:
            self.logger.error(f"Failed {self.operation} after {duration:.2f}s: {exc_val}")
            # Don't suppress the exception
            return False

if __name__ == "__main__":
    # Test the error handling utilities
    logger = setup_logging()
    
    @error_handler
    def test_function():
        """Test function to demonstrate error handling."""
        validate_input_files("/tmp/test_file.txt")
        validate_output_directory("/tmp/test_output")
        
        # Test safe operations
        result = safe_division(10, 0, default=1.0)
        percentage = safe_percentage(50, 100)
        
        logger.info(f"Safe division result: {result}")
        logger.info(f"Percentage: {percentage}%")
        
        # Test error context
        with ErrorContext("test operation", logger):
            # Simulate some work
            import time
            time.sleep(0.1)
        
        return "Success"
    
    try:
        result = test_function()
        logger.info(f"Test completed: {result}")
    except Exception as e:
        logger.error(f"Test failed: {e}")