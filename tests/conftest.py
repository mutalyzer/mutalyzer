import socket


class BlockSocket(socket.socket):
    def __init__(self, *args, **kwargs):
        raise Exception("Socket blocked!")


socket.socket = BlockSocket
