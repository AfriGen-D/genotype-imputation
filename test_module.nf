process TEST_PROCESS {
    tag "test"
    label 'process_single'
    
    script:
    """
    echo "test"
    """
}
