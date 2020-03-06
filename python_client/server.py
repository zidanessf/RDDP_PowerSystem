import testsvc_pb2,testsvc_pb2_grpc,grpc
from concurrent import futures
import time
class Dispatcher(testsvc_pb2_grpc.TestMathServicer):
	def Mul(self, request, context):
		# missing associated documentation comment in .proto file
		return testsvc_pb2.BinaryOpResp(result=1)

	def Add(self, request, context):
		# missing associated documentation comment in .proto file
		print(1)
		return testsvc_pb2.BinaryOpResp(result=1)
def serve():
	server = grpc.server(futures.ThreadPoolExecutor(max_workers=10))
	testsvc_pb2_grpc.add_TestMathServicer_to_server(Dispatcher(), server)
	server.add_insecure_port('127.0.0.1:19999')
	server.start()
	print('server started in port 19999')
	try:
		while True:
			time.sleep(1)
	except KeyboardInterrupt:
		server.stop(0)

if __name__ == '__main__':
	serve()