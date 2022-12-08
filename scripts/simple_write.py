def simple_write(filename, content=None):
    '''
    Writes a string to a file, handling opening and closing.
    If parameter conent is not supplied, acts as a touch.
    This function overwrites.

    :param filename: str, the path to the file we are creating and writing.
    :param content: str (optional), what we are writing
    '''
    with open(filename, 'w') as f:
        if content is not None:
            f.write(content)