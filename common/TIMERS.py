import time

class Timer:
    message_start = 'Starting task.'
    message_finished = 'Finished! Yay!.'
    message_elapsed_time = "Time elapsed: "

    def __init__(self, message_start=None, message_finished=None, message_elapsed=None) -> None:
        if message_start is not None:
            self.message_start = message_start
        if message_finished is not None:
            self.message_finished = message_finished
        if message_elapsed is not None:
            self.message_elapsed = message_elapsed

    def _get_duration(self, function, *args, **kwargs) -> float:
        t1 = time.time()
        function(args, kwargs)
        return time.time() - t1

    def print_duration(self, function, *args, **kwargs):
        self.print_message(self.message_start)      
        time_total = self._get_duration(function, args, kwargs)
        self.print_message(self.message_elapsed_time + str(time_total))
        self.print_message(self.message_finished)

    def print_elapsed_time_message(self, time):
        print(f"Time elapsed: {time}")

    def print_message(self, message):
        print(message)

